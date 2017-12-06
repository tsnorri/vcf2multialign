/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <vcf2multialign/channel_sink.hh>
#include <vcf2multialign/dispatch_fn.hh>


// FIXME: at the moment, this file needs to be compiled with Clang.
namespace vcf2multialign {
	
	void channel_sink_impl::open(int fd, dispatch_ptr <dispatch_semaphore_t> const &writing_semaphore)
	{
		auto const qname(boost::str(boost::format("fi.iki.tsnorri.vcf2multialign.channel-queue-%u") % fd));
		dispatch_ptr <dispatch_queue_t> channel_queue(dispatch_queue_create(qname.c_str(), DISPATCH_QUEUE_SERIAL));
		dispatch_ptr <dispatch_io_t> channel(dispatch_io_create(DISPATCH_IO_STREAM, fd, *channel_queue, ^(int const error) {
			if (error)
				std::cerr << "IO channel for file " << fd << " closed unexpectedly: " << error << std::endl;
		}));
		
		m_channel_queue = std::move(channel_queue);
		m_channel = std::move(channel);
		m_writing_semaphore = writing_semaphore;
		m_writing_group.reset(dispatch_group_create());
	}
	
	
	void channel_sink_impl::close()
	{
		auto const close_flag(m_close_flag.exchange(true));
		always_assert(!close_flag);
		
		// Check for closing group.
		auto closing_group(*m_closing_group);
		if (closing_group)
			dispatch_group_enter(closing_group);
		
		// Wait for the writing operations to be ready.
		dispatch_group_notify(*m_writing_group, *m_channel_queue, ^{
			dispatch_io_close(*m_channel, 0);
			
			if (closing_group)
				dispatch_group_leave(closing_group);
			
			// All instance variables are valid up to this point.
			delete(this);
		});
	}
	
	
	std::streamsize channel_sink_impl::write(char const *bytes, std::streamsize const n)
	{
		auto buffer(dispatch_data_create(bytes, n, *m_channel_queue, DISPATCH_DATA_DESTRUCTOR_DEFAULT));
		dispatch_semaphore_wait(*m_writing_semaphore, DISPATCH_TIME_FOREVER);
		
		dispatch_group_enter(*m_writing_group);
		// Check that the sink has not been closed while waiting.
		// (This is still unsafe b.c. the function calls may be interleaved s.t. m_close_flag is set
		// just after the check below. The intent is to find some clear bugs. This could be made safe
		// by replacing m_close_flag with an ordinary bool, adding a dispatch_semaphore_t with value one
		// and waiting on the semaphore both here and in close().)
		always_assert(!m_close_flag);
		dispatch_io_write(*m_channel, 0, buffer, *m_channel_queue, ^(bool const done, dispatch_data_t, int const error) {
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
		
		return n;
	}
}
