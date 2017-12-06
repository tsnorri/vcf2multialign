/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <vcf2multialign/gzip_sink.hh>
#include <vcf2multialign/dispatch_fn.hh>


// FIXME: at the moment, this file needs to be compiled with Clang.
namespace vcf2multialign {
	
	gzip_sink_impl::pool_resource_type *gzip_sink_impl::s_pool_resource = nullptr;
	
	
	unsigned char *gzip_sink_impl::allocate_input_buffer(std::size_t const n)
	{
		return reinterpret_cast <unsigned char *>(s_pool_resource->allocate(n, alignof(unsigned char)));
	}
	
	
	void gzip_sink_impl::deallocate_input_buffer(unsigned char *buffer, std::size_t const n)
	{
		s_pool_resource->deallocate(buffer, n, alignof(unsigned char));
	}
	
	
	unsigned char *gzip_sink_impl::allocate_compression_buffer()
	{
		return reinterpret_cast <unsigned char *>(s_pool_resource->allocate(m_write_size, alignof(unsigned char)));
	}
	
	
	void gzip_sink_impl::deallocate_compression_buffer(unsigned char *buffer)
	{
		s_pool_resource->deallocate(buffer, m_write_size, alignof(unsigned char));
	}
	
	
	void gzip_sink_impl::allocate_and_place_compression_buffer()
	{
		auto *dst(allocate_compression_buffer());
		m_compression_stream.avail_out	= m_write_size;
		m_compression_stream.next_out	= dst;
	}

	
	void gzip_sink_impl::close()
	{
		auto const close_flag(m_close_flag.exchange(true));
		always_assert(!close_flag);
		
		// Check for closing group.
		auto closing_group(*m_closing_group);
		if (closing_group)
			dispatch_group_enter(closing_group);
		
		// Place an operation after pending compression operations.
		dispatch_async(*m_compression_queue, ^{
			// Start writing operations for the remaining compressed data.
			while (true)
			{
				auto const st(deflate(&m_compression_stream, Z_FINISH));
				if (Z_STREAM_END != st)
					write_compressed_data_and_replace_buffer(true);
				else
				{
					write_compressed_data_and_replace_buffer(false);
					break;
				}
			}
			
			// Free allocated memory.
			deflateEnd(&m_compression_stream);
			
			// Wait for the writing operations to be ready.
			dispatch_group_notify(*m_writing_group, *m_channel_queue, ^{
				dispatch_io_close(*m_channel, 0);
			
				if (closing_group)
					dispatch_group_leave(closing_group);
			
				// All instance variables are valid up to this point.
				delete(this);
			});
		});
	}
	
	
	void gzip_sink_impl::open(
		int fd,
		dispatch_ptr <dispatch_semaphore_t> const &compression_semaphore,
		dispatch_ptr <dispatch_semaphore_t> const &writing_semaphore
	)
	{
		// Create the queues and IO channel.
		auto const compress_qname(boost::str(boost::format("fi.iki.tsnorri.vcf2multialign.gzip-channel-compress-queue-%u") % fd));
		auto const write_qname(boost::str(boost::format("fi.iki.tsnorri.vcf2multialign.gzip-channel-write-queue-%u") % fd));
		dispatch_ptr <dispatch_queue_t> compression_queue(dispatch_queue_create(compress_qname.c_str(), DISPATCH_QUEUE_SERIAL));
		dispatch_ptr <dispatch_queue_t> channel_queue(dispatch_queue_create(write_qname.c_str(), DISPATCH_QUEUE_SERIAL));
		dispatch_ptr <dispatch_io_t> channel(dispatch_io_create(DISPATCH_IO_STREAM, fd, *channel_queue, ^(int const error) {
			if (error)
			{
				std::cerr << "IO channel for file " << fd << " closed unexpectedly: " << error << std::endl;
				throw std::runtime_error("");
			}
		}));
		
		m_write_size = 65536; // FIXME: come up with a better guess.
		
		// Initialize variables.
		m_compression_queue = std::move(compression_queue);
		m_channel_queue = std::move(channel_queue);
		m_channel = std::move(channel);
		m_compression_semaphore = compression_semaphore;
		m_writing_semaphore = writing_semaphore;
		m_writing_group.reset(dispatch_group_create());
		
		// Initialize m_compression_stream.
		{
			m_compression_stream.zalloc = Z_NULL;
			m_compression_stream.zfree = Z_NULL;
			auto const st(deflateInit2(
				&m_compression_stream,
				Z_DEFAULT_COMPRESSION,
				Z_DEFLATED,
				16 + 15,
				8,
				Z_DEFAULT_STRATEGY
			));
			always_assert(Z_OK == st);
			allocate_and_place_compression_buffer();
		}
	}
	
	
	void gzip_sink_impl::write_compressed_data_and_replace_buffer(bool get_new_buffer)
	{
		dispatch_semaphore_wait(*m_writing_semaphore, DISPATCH_TIME_FOREVER);
		
		// Take ownership of the compression buffer.
		auto cs_buffer(m_compression_stream.next_out - (m_write_size - m_compression_stream.avail_out));
		auto dispatch_buffer(dispatch_data_create(
			cs_buffer,
			m_write_size - m_compression_stream.avail_out,
			*m_channel_queue,
			^{ deallocate_compression_buffer(cs_buffer); }
		));
		
		// Get a new buffer for deflate.
		if (get_new_buffer)
			allocate_and_place_compression_buffer();
		else
		{
			m_compression_stream.next_out = nullptr;
			m_compression_stream.avail_out = 0;
		}
		
		// Write the contents of the dispatch buffer.
		dispatch_group_enter(*m_writing_group);
		dispatch_io_write(*m_channel, 0, dispatch_buffer, *m_channel_queue, ^(bool const done, dispatch_data_t, int const error) {
			if (error)
			{
				auto const fd(dispatch_io_get_descriptor(*m_channel));
				std::cerr << "IO channel failed to write to file " << fd << '.' << std::endl;
				abort();
			}
			
			if (done)
			{
				dispatch_semaphore_signal(*m_writing_semaphore);
				dispatch_group_leave(*m_writing_group);
			}
		});
	}
	
	
	std::streamsize gzip_sink_impl::write(char const *bytes, std::streamsize const n)
	{
		// Wait on the sema.
		dispatch_semaphore_wait(*m_compression_semaphore, DISPATCH_TIME_FOREVER);
		
		dispatch_group_enter(*m_writing_group);
		always_assert(!m_close_flag);
		
		// Copy the data for compression_queue. The block below will own input.
		auto *input(allocate_input_buffer(n));
		memcpy(input, bytes, n);
		dispatch_async(*m_compression_queue, ^{
			// Compress the data and write in blocks.
			m_compression_stream.avail_in	= n;
			m_compression_stream.next_in	= input;
			
			while (m_compression_stream.avail_in)
			{
				assert(m_compression_stream.next_out);
				
				// Compress the input.
				auto const st(deflate(&m_compression_stream, Z_NO_FLUSH));
				always_assert(Z_OK == st);
				
				// Check the output buffer.
				if (0 == m_compression_stream.avail_out)
				{
					// Create a dispatch data object from the current buffer.
					write_compressed_data_and_replace_buffer(true);
				}
			}
			
			deallocate_input_buffer(input, n);
			dispatch_semaphore_signal(*m_compression_semaphore);
			dispatch_group_leave(*m_writing_group);
		});
		
		return n;
	}
}
