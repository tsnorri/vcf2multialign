/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_BIT_STRING_FILE_HH
#define VCF2MULTIALIGN_BIT_STRING_FILE_HH

#include <vcf2multialign/types.hh>


namespace vcf2multialign {

	class bit_string_file
	{
	protected:
		file_ostream	*m_stream{};
	
	public:
		bit_string_file() = default;

		bit_string_file(file_ostream &stream):
			m_stream(&stream)
		{
		}

		// Just write bytes for now.
		bool is_open() const				{ return m_stream && m_stream->is_open(); }
		void close()						{ if (m_stream) m_stream->close(); m_stream = nullptr; }
		void write_ones(std::size_t count)	{ std::fill_n(std::ostream_iterator <char>(*m_stream), count, '1'); }
		void write_zeros(std::size_t count)	{ std::fill_n(std::ostream_iterator <char>(*m_stream), count, '0'); }
		void flush()						{ *m_stream << std::flush; }
	};
}

#endif
