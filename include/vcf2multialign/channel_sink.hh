/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CHANNEL_SINK_HH
#define VCF2MULTIALIGN_CHANNEL_SINK_HH

#include <boost/iostreams/categories.hpp>
#include <iosfwd>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/polymorphic_sink.hh>
#include <vcf2multialign/util.hh>


namespace vcf2multialign {
	
	class channel_sink_impl final : public polymorphic_sink_impl
	{
	protected:
		dispatch_ptr <dispatch_io_t>		m_channel{};
		dispatch_ptr <dispatch_queue_t>		m_channel_queue{};
		dispatch_ptr <dispatch_semaphore_t>	m_writing_semaphore{};
		dispatch_ptr <dispatch_group_t>		m_writing_group{};
		dispatch_ptr <dispatch_group_t>		m_closing_group{};
		std::atomic_bool					m_close_flag{false};
		
	public:
		void close();
		void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) { m_closing_group = group; }
		void open(int fd, dispatch_ptr <dispatch_semaphore_t> const &writing_semaphore);
		std::streamsize write(char const *bytes, std::streamsize const n);
	};
	
	
	class channel_sink final
	{
	public:
		typedef char						char_type;
		struct category :
			boost::iostreams::sink_tag,
			boost::iostreams::closable_tag
		{
		};
		
	protected:
		// FIXME: copying should work with our current implementation but come up with a better way.
		channel_sink_impl	*m_impl{};
		
	public:
		channel_sink() = default;
		
		void close() { if (m_impl) m_impl->close(); }
		void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) { m_impl->set_closing_group(group); }
		void open(int fd, dispatch_ptr <dispatch_semaphore_t> const &writing_semaphore) { m_impl = new channel_sink_impl; m_impl->open(fd, writing_semaphore); }
		std::streamsize write(char const *bytes, std::streamsize const n) { return m_impl->write(bytes, n); }
	};
}
	
#endif
