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
		
		m_d.channel_queue = std::move(channel_queue);
		m_d.channel = std::move(channel);
		m_d.write_semaphore = std::move(write_semaphore);
	}
	
	
	void channel_sink::close()
	{
		auto const close_flag(m_close_flag.exchange(true));
		always_assert(!close_flag);

		// set_closing_group needs to be called in the same queue as close().
		dispatch_group_t closing_group(*m_d.closing_group);
		dispatch_queue_t channel_queue(*m_d.channel_queue);
		auto fn = [d = std::move(m_d)]() mutable {
			dispatch_io_close(*d.channel, 0);
		};
		
		// fn now owns m_d contents.
		// Make sure that all the queued write operations have been completed before closing
		// and possibly notify the dispatch group.
		if (closing_group)
			dispatch_group_async_fn(closing_group, channel_queue, std::move(fn));
		else
			dispatch_async_fn(channel_queue, std::move(fn));
	}
	
	
	std::streamsize channel_sink::write(char const *bytes, std::streamsize const n)
	{
		auto queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
		auto buffer(dispatch_data_create(bytes, n, queue, DISPATCH_DATA_DESTRUCTOR_DEFAULT));
		
		dispatch_semaphore_t semaphore(*m_d.write_semaphore);
		dispatch_retain(semaphore);
		int const fd(dispatch_io_get_descriptor(*this->m_d.channel));
		// In case other variables are needed in the block, retain them here.
		dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
		
		// Check that the sink has not been closed while waiting.
		// (This is still unsafe b.c. the function calls may be interleaved s.t. m_close_flag is set
		// just after the check below. The intent is to find some clear bugs. This could be made safe
		// by replacing m_close_flag with an ordinary bool, adding a dispatch_semaphore_t with value one
		// and waiting on the semaphore both here and in close().)
		always_assert(!m_close_flag);
		dispatch_io_write(*m_d.channel, 0, buffer, *m_d.channel_queue, ^(bool const done, dispatch_data_t, int const error) {
			if (error)
			{
				std::cerr << "IO channel failed to write to file " << fd << '.' << std::endl;
				abort();
			}
			
			if (done)
			{
				dispatch_semaphore_signal(semaphore);
				dispatch_release(semaphore);
			}
		});
		
		return n;
	}
}
