/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCF_READER_HH
#define VCF2MULTIALIGN_VCF_READER_HH

#include <istream>
#include <map>
#include <set>
#include <experimental/string_view>
#include <stdexcept>
#include <vector>


// XXX Hack.
namespace std {
	using std::experimental::string_view;
}


namespace  vcf2multialign { namespace detail {
	template <typename t_enum>
	constexpr typename std::underlying_type <t_enum>::type to_underlying(t_enum e)
	{
		return static_cast <typename std::underlying_type <t_enum>::type>(e);
	}

	
	template <typename t_integral>
	void parse_int(std::string_view const &str, t_integral &res)
	{
		res = 0;
		for (auto const c : str)
		{
			if (! ('0' <= c && c <= '9'))
				throw std::runtime_error("Unable to convert string into number");
			
			res *= 10;
			res += c - '0';
		}
	}
	
	
	// Split a t_string into string_views.
	template <bool t_emplace_back, typename t_string, typename t_invalid_pos>
	size_t read_fields(
		t_string const &string,
		t_invalid_pos const invalid_pos,
		char const *string_start,
		char const *sep,
		size_t const count,
		std::vector <std::string_view> &res
	)
	{
		char const *string_ptr(string_start);
		
		if (t_emplace_back)
			res.clear();
		else
			res.resize(count);
		
		typename t_string::size_type sep_pos(0);
		size_t i(0);
		bool should_break(false);
		while (i < count)
		{
			sep_pos = string.find(sep, sep_pos);
			if (invalid_pos == sep_pos)
			{
				should_break = true;
				sep_pos = string.size();
			}
			
			if (t_emplace_back)
				res.emplace_back(string_ptr, sep_pos - (string_ptr - string_start));
			else
			{
				std::string_view sv(string_ptr, sep_pos - (string_ptr - string_start));
				res[i] = std::move(sv);
			}
			
			++i;
			++sep_pos;
			string_ptr = string_start + sep_pos;
			
			if (should_break)
				break;
		}
		return i;
	}
							   
	
	// Split a string into string_views.
	template <bool t_emplace_back>
	size_t read_fields(
		std::string const &string,
		char const *sep,
		size_t const count,
		std::vector <std::string_view> &res
	)
	{
		return read_fields <t_emplace_back>(string, std::string::npos, string.c_str(), sep, count, res);
	}
	
	template <bool t_emplace_back>
	size_t read_fields(
		std::string_view const &sv,
		char const *sep,
		size_t const count,
		std::vector <std::string_view> &res
	)
	{
		return read_fields <t_emplace_back>(sv, std::string_view::npos, static_cast <char const *>(sv.data()), sep, count, res);
	}
}}


namespace vcf2multialign {
	
	enum class vcf_field : uint8_t {
		CHROM	= 0,
		POS		= 1,
		ID		= 2,
		REF		= 3,
		ALT		= 4,
		QUAL	= 5,
		FILTER	= 6,
		INFO	= 7,
		FORMAT	= 8,
		ALL		= 9
	};
	
	
	enum { NULL_ALLELE = std::numeric_limits <uint8_t>::max() };
	
	
	class vcf_reader;
	
	
	class variant
	{
		friend class vcf_reader;
		
	public:
		typedef std::vector <std::string_view>	sample_field_vector;
		typedef std::vector <uint8_t>			genotype_vector;
		
	protected:
		std::vector <std::string_view> m_var_fields;
		std::vector <std::string_view> m_alt;
		std::string m_format;
		std::map <std::string, uint8_t> m_format_fields;
		std::set <std::string> m_requested_format_fields;
		std::size_t m_pos{0};
		std::size_t m_sample_count{0};
		std::size_t m_lineno{0};
		uint8_t m_format_max{0};
		bool m_parsed_alt{false};
		
	public:
		variant(std::size_t const sample_count, std::set <std::string> const &requested_format_fields):
			m_var_fields(sample_count),
			m_requested_format_fields(requested_format_fields),
			m_sample_count(sample_count)
		{
		}
		
		variant(std::size_t const sample_count = 0):
			m_var_fields(sample_count),
			m_sample_count(sample_count)
		{
		}
		
		void add_format_field(std::string const &name) { m_requested_format_fields.emplace(name); }
		
		std::size_t const lineno() const		{ return m_lineno; }
		std::string_view const &chrom() const	{ return m_var_fields[detail::to_underlying(vcf_field::CHROM)]; }
		std::string_view const &id() const		{ return m_var_fields[detail::to_underlying(vcf_field::ID)]; }
		std::string_view const &ref() const		{ return m_var_fields[detail::to_underlying(vcf_field::REF)]; }
		
		size_t pos();
		size_t zero_based_pos();
		std::vector <std::string_view> const &alt();
		void reset();
		
		// Check the format of the current sample.
		void prepare_samples();
		
		// Split the sample into fields up to the last one in requested_format_fields, thread-safe.
		void parse_sample(size_t const sample_no, sample_field_vector /* out */ &sample_fields) const;
		
		// Get the genotype from a parsed sample.
		void get_genotype(sample_field_vector const &sample_fields, genotype_vector /* out */ &res, bool /* out */ &phased) const;
		
	protected:
		void map_format_fields(std::string_view const &format);
	};
	
	
	class vcf_reader
	{
	public:
		typedef std::map <std::string, std::size_t> sample_name_map;
		
	protected:
		std::istream *m_stream{nullptr};
		sample_name_map m_sample_names;
		std::string m_line;
		std::size_t m_lineno{0};
		size_t m_parsed_field_count{0};
		
		std::istream::pos_type m_first_variant_offset{0};
		std::size_t m_last_header_lineno{0};
		
	public:
		vcf_reader(std::istream &stream):
			m_stream(&stream)
		{
		}
		
		vcf_reader() {}
		
		void set_stream(std::istream &stream) { m_stream = &stream; }
		std::size_t lineno() const { return m_lineno; }
		void read_header();
		void reset();
		bool get_line(std::string &line);
		void get_next_variant(variant &var, std::string &line) const;
		bool get_next_variant(variant &var);
		size_t sample_no(std::string const &sample_name) const;
		size_t sample_count() const { return m_sample_names.size(); }
		void set_parsed_fields(vcf_field const last_field);
		sample_name_map const &sample_names() { return m_sample_names; }
	};
}

#endif
