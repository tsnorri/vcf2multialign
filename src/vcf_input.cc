/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/vcf_input.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace ios = boost::iostreams;


namespace vcf2multialign {
	
	void vcf_stream_input_base::reset_to_first_variant_offset()
	{
		stream_reset();
		m_len = 0;
		m_pos = 0;
	}
	
	
	bool vcf_stream_input_base::getline(std::string_view &dst)
	{
		bool const st(stream_getline());
		if (!st)
			return false;
		
		dst = std::string_view(m_buffer.data(), m_buffer.size());
		return true;
	}
	
	
	void vcf_stream_input_base::store_first_variant_offset()
	{
		m_first_variant_offset = stream_tellg();
	}
	
	
	void vcf_stream_input_base::fill_buffer(vcf_reader &reader)
	{
		// Copy the remainder to the beginning.
		if (m_pos + 1 < m_len)
		{
			char *data_start((char *) m_buffer.data()); // FIXME: hack, my libc++ does not yet have the non-const data().
			char *start(data_start + m_pos + 1);
			char *end(data_start + m_len);
			std::copy(start, end, data_start);
			m_len -= m_pos + 1;
		}
		else
		{
			m_len = 0;
		}
		
		// Read until there's at least one newline in the buffer.
		while (true)
		{
			char *data_start((char *) m_buffer.data()); // FIXME: hack, my libc++ does not yet have the non-const data().
			char *data(data_start + m_len);
		
			std::size_t space(m_buffer.size() - m_len);
			stream_read(data, space);
			std::size_t const read_len(stream_gcount());
			m_len += read_len;
		
			if (stream_eof())
			{
				m_pos = m_len;
				auto const end(data_start + m_len);
				reader.set_p(data_start);
				reader.set_pe(end);
				reader.set_eof(end);
				return;
			}
		
			// Try to find the last newline in the new part.
			std::string_view sv(data, read_len);
			m_pos = sv.rfind('\n');
			if (std::string_view::npos != m_pos)
			{
				m_pos += (data - data_start);
				reader.set_p(data_start);
				reader.set_pe(data_start + m_pos + 1);
				return;
			}
		
			m_buffer.resize(2 * m_buffer.size());
		}
	}
	
	
	bool vcf_mmap_input::getline(std::string_view &dst)
	{
		auto const size(m_handle.size());
		assert(m_pos <= size);
		
		auto const sv(m_handle.operator std::string_view());
		auto const nl_pos(sv.find('\n', m_pos));
		if (std::string_view::npos == nl_pos)
			return false;
		
		// Found the newline, update the variables.
		dst = sv.substr(m_pos, nl_pos - m_pos);
		m_pos = 1 + nl_pos;
		return true;
	}
	
	
	void vcf_mmap_input::store_first_variant_offset()
	{
		m_first_variant_offset = m_pos;
	}
	
	
	void vcf_mmap_input::fill_buffer(vcf_reader &reader)
	{
		auto const *bytes(m_handle.data());
		auto const len(m_handle.size());
		reader.set_p(bytes + m_first_variant_offset);
		reader.set_pe(bytes + len);
		reader.set_eof(bytes + len);
	}
}
