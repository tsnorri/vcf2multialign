/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <fcntl.h>
#include <iostream>
#include <sys/stat.h>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {

	size_t variant::pos()
	{
		if (0 == m_pos)
			detail::parse_int(m_var_fields[detail::to_underlying(vcf_field::POS)], m_pos);
		
		return m_pos;
	}
	
	
	size_t variant::zero_based_pos()
	{
		auto const p(pos());
		if (0 == p)
			throw std::runtime_error("Unexpected position");
		
		return p - 1;
	}
	
	
	std::vector <std::string_view> const &variant::alt()
	{
		if (!m_parsed_alt)
			detail::read_fields <true>(m_var_fields[detail::to_underlying(vcf_field::ALT)], ",", -1, m_alt);
		
		return m_alt;
	}
	
	
	void variant::reset()
	{
		m_alt.clear();
		m_pos = 0;
		m_parsed_alt = false;
		m_format = "";
		m_format_fields.clear();
		m_format_max = 0;
		m_lineno = 0;
	}
	
	
	void variant::map_format_fields(std::string_view const &format)
	{
		m_format_fields.clear();
		auto format_fields(m_requested_format_fields);
		char const *string_start(static_cast <char const *>(format.data()));
		char const *string_ptr(string_start);
		uint8_t i(0);
		
		std::string_view::size_type sep_pos(0);
		std::string buffer;
		bool should_break(false);
		while (true)
		{
			sep_pos = format.find(":", sep_pos);
			if (std::string_view::npos == sep_pos)
			{
				should_break = true;
				sep_pos = format.size();
			}
			std::string_view sv(string_ptr, sep_pos - (string_ptr - string_start));
			
			buffer.assign(sv.cbegin(), sv.cend());
			auto it(format_fields.find(buffer));
			if (format_fields.cend() != it)
			{
				m_format_fields.insert(std::make_pair(buffer, i));
				m_format_max = i;
			}
			
			++i;

			if (0 == format_fields.size())
				break;
			
			if (should_break)
				break;
			
			string_ptr = string_start + sep_pos;
			++sep_pos;
		}
		
		m_format.assign(format.cbegin(), format.cend());
	}
	
	
	// Parse the format field if needed.
	void variant::prepare_samples()
	{
		auto const &format(m_var_fields[detail::to_underlying(vcf_field::FORMAT)]);
		if (m_format != format)
			map_format_fields(format);
	}

	
	// Parse a sample according to the previously read format.
	void variant::parse_sample(size_t const sample_no, std::vector <std::string_view> /* out */ &sample_fields) const
	{
		auto const idx(detail::to_underlying(vcf_field::ALL) + sample_no - 1);
		auto const &val(m_var_fields[idx]);
		detail::read_fields <false>(val, ":", 1 + m_format_max, sample_fields);
	}

	
	// Read the genotype field from a parsed sample.
	void variant::get_genotype(std::vector <std::string_view> const &sample_fields, std::vector <uint8_t> /* out */ &res, bool /* out */ &phased) const
	{
		res.clear();
		uint8_t alt_num(0);
		phased = true;
		
		auto const gt_idx_it(m_format_fields.find("GT"));
		if (m_format_fields.cend() == gt_idx_it)
			throw std::runtime_error("GT not present in format");
		
		bool expect_separator(false);
		for (auto c : sample_fields[gt_idx_it->second])
		{
			if ('.' == c)
			{
				alt_num = NULL_ALLELE;

				if (expect_separator)
					throw std::runtime_error("Expected separator");
				
				expect_separator = true;
			}
			else if ('0' <= c && c <= '9')
			{
				alt_num *= 10;
				alt_num += c - '0';

				if (expect_separator)
					throw std::runtime_error("Expected separator");
			}
			else
			{
				if ('/' == c)
					phased = false;
				else if ('|' != c)
					throw std::runtime_error("Unexpected character");
				
				res.emplace_back(alt_num);
				alt_num = 0;
				expect_separator = false;
			}
		}
		
		res.emplace_back(alt_num);
	}
	
	
	// Parse the VCF header.
	void vcf_reader::read_header()
	{
		std::vector <std::string_view> fields(m_parsed_field_count);
		
		// For now, just skip lines that begin with "##".
		while (std::getline(*m_stream, m_line))
		{
			++m_lineno;
			if (! ('#' == m_line[0] && '#' == m_line[1]))
				break;
			
			// FIXME: header handling goes here.
		}
		
		// Check for column headers.
		{
			auto const prefix(std::string("#CHROM"));
			auto const res(mismatch(m_line.cbegin(), m_line.cend(), prefix.cbegin(), prefix.cend()));
			if (prefix.cend() != res.second)
			{
				std::cerr << "Expected the next line to start with '#CHROM', instead got: '" << m_line << "'" << std::endl;
				exit(1);
			}
		}
		
		// Read sample names.
		std::vector <std::string_view> header_names;
		detail::read_fields <true>(m_line, "\t", -1, header_names);
		std::size_t i(1);
		for (auto it(header_names.cbegin() + detail::to_underlying(vcf_field::ALL)), end(header_names.cend()); it != end; ++it)
		{
			std::string str(*it);
			auto const res(m_sample_names.emplace(std::move(str), i));
			if (!res.second)
				throw std::runtime_error("Duplicate sample name");
			++i;
		}
		
		// stream now points to the first variant.
		m_first_variant_offset = m_stream->tellg();
		m_last_header_lineno = m_lineno;
	}
	
	
	// Seek to the beginning of the records.
	void vcf_reader::reset()
	{
		m_stream->clear();
		m_stream->seekg(m_first_variant_offset);
		m_lineno = m_last_header_lineno;
	}
	
	
	// Fill var with the first m_parsed_field_count fields.
	bool vcf_reader::get_next_variant(variant &var)
	{
		if (!std::getline(*m_stream, m_line))
			return false;
		
		++m_lineno;
		var.reset();
		detail::read_fields <false>(m_line, "\t", m_parsed_field_count, var.m_var_fields);
		var.m_lineno = m_lineno;
		return true;
	}
	
	
	// Return the 1-based number of the given sample.
	size_t vcf_reader::sample_no(std::string const &sample_name) const
	{
		auto const it(m_sample_names.find(sample_name));
		if (it == m_sample_names.cend())
			return 0;
		return it->second;
	}
	
	
	// Set the last field to be parsed on a line.
	void vcf_reader::set_parsed_fields(vcf_field const last_field)
	{
		if (last_field == vcf_field::ALL)
			m_parsed_field_count = sample_count() + detail::to_underlying(vcf_field::ALL);
		else
			m_parsed_field_count = 1 + detail::to_underlying(last_field);
	}
}
