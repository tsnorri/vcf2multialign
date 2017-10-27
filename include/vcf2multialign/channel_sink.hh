/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CHANNEL_SINK_HH
#define VCF2MULTIALIGN_CHANNEL_SINK_HH

#include <boost/iostreams/categories.hpp>
#include <iosfwd>
#include <vcf2multialign/dispatch_fn.hh>


namespace vcf2multialign {
	
	class channel_sink
	{
	public:
		typedef char						char_type;
		typedef boost::iostreams::sink_tag	category;
		
	protected:
		dispatch_ptr <dispatch_io_t>		m_channel{};
		dispatch_ptr <dispatch_queue_t>		m_channel_queue{};
		dispatch_ptr <dispatch_semaphore_t>	m_write_semaphore{};
		
	public:
		channel_sink() = default;
		
		void close();
		void close_async(dispatch_group_t group);
		void open(int fd, dispatch_ptr <dispatch_semaphore_t> const &write_semaphore);
		std::streamsize write(char const *bytes, std::streamsize const n);
		
	protected:
		void close_2();
	};
}
	
#endif
