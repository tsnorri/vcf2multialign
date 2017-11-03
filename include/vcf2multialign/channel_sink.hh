/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CHANNEL_SINK_HH
#define VCF2MULTIALIGN_CHANNEL_SINK_HH

#include <boost/iostreams/categories.hpp>
#include <iosfwd>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/util.hh>


namespace vcf2multialign {
	
	class channel_sink
	{
	public:
		typedef char						char_type;
		typedef boost::iostreams::sink_tag	category;
		
	protected:
		// Use a container to store the instance variables when closing.
		struct data
		{
			dispatch_ptr <dispatch_io_t>		channel{};
			dispatch_ptr <dispatch_queue_t>		channel_queue{};
			dispatch_ptr <dispatch_semaphore_t>	write_semaphore{};
			dispatch_ptr <dispatch_group_t>		closing_group{};
		};
		
	protected:
		data					m_d;
		copyable_atomic <bool>	m_close_flag{false};
		
	public:
		channel_sink() = default;
		
		void close();
		void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) { m_d.closing_group = group; }
		void open(int fd, dispatch_ptr <dispatch_semaphore_t> const &write_semaphore);
		std::streamsize write(char const *bytes, std::streamsize const n);
		
	protected:
		void close_2();
	};
}
	
#endif
