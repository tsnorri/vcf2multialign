/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <vcf2multialign/channel_sink.hh>
#include <vcf2multialign/dispatch_fn.hh>


// FIXME: at the moment, this file needs to be compiled with Clang.
namespace vcf2multialign {
	
	void channel_sink::open(int fd, dispatch_ptr <dispatch_semaphore_t> const &write_semaphore)
	{
		auto const qname(boost::str(boost::format("fi.iki.tsnorri.vcf2multialign.channel-queue-%u") % fd));
		dispatch_ptr <dispatch_queue_t> channel_queue(dispatch_queue_create(qname.c_str(), DISPATCH_QUEUE_SERIAL));
		dispatch_ptr <dispatch_io_t> channel(dispatch_io_create(DISPATCH_IO_STREAM, fd, *channel_queue, ^(int const error) {
			if (error)
				std::cerr << "IO channel for file " << fd << " closed unexpectedly: " << error << std::endl;
		}));
		
		m_channel_queue = std::move(channel_queue);
		m_channel = std::move(channel);
		m_write_semaphore = std::move(write_semaphore);
	}
	
	
	void channel_sink::close_async(dispatch_group_t group)
	{
		// Make sure that all the queued write operations have been completed before closing
		// and notify the dispatch group.
		dispatch(this).group_async <&channel_sink::close_2>(group, *m_channel_queue);
	}
	
	
	void channel_sink::close()
	{
		// Make sure that all the queued write operations have been completed before closing
		// and wait in the current thread.
		dispatch(this).sync <&channel_sink::close_2>(*m_channel_queue);
	}
	
	
	void channel_sink::close_2()
	{
		// Safe b.c. close_2 is called in m_channel_queue.
		auto channel(*m_channel);
		if (channel)
		{
			dispatch_io_close(channel, 0);
			m_channel.reset();
		}
	}
	
	
	std::streamsize channel_sink::write(char const *bytes, std::streamsize const n)
	{
		auto queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
		auto buffer(dispatch_data_create(bytes, n, queue, DISPATCH_DATA_DESTRUCTOR_DEFAULT));
		
		dispatch_semaphore_wait(*m_write_semaphore, DISPATCH_TIME_FOREVER);
		dispatch_io_write(*m_channel, 0, buffer, *m_channel_queue, ^(bool const done, dispatch_data_t, int error) {
			if (error)
			{
				std::cerr << "IO channel failed to write to file " << dispatch_io_get_descriptor(*this->m_channel) << '.' << std::endl;
				abort();
			}
			
			if (done)
				dispatch_semaphore_signal(*this->m_write_semaphore);
		});
		
		return n;
	}
}
