/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_POLYMORPHIC_SINK_HH
#define VCF2MULTIALIGN_POLYMORPHIC_SINK_HH

#include <boost/iostreams/categories.hpp>


namespace vcf2multialign {
	
	struct polymorphic_sink_impl
	{
		virtual ~polymorphic_sink_impl() {}
		virtual void close() = 0;
		virtual void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) {}
		virtual std::streamsize write(char const *bytes, std::streamsize const n) = 0;
	};
	
	
	class polymorphic_sink
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
		polymorphic_sink_impl	*m_impl{};
		
	public:
		void open(polymorphic_sink_impl *impl) { m_impl = impl; }
		void close() { if (m_impl) m_impl->close(); }
		void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) { assert(m_impl); m_impl->set_closing_group(group); }
		std::streamsize write(char const *bytes, std::streamsize const n) { assert(m_impl); return m_impl->write(bytes, n); }
	};
}

#endif
