/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCF_READER_HH
#define VCF2MULTIALIGN_VCF_READER_HH

#include <istream>
#include <map>
#include <vector>
#include <vcf2multialign/variant.hh>


namespace vcf2multialign {
	
	class vcf_reader
	{
	public:
		typedef std::function <bool(transient_variant const &var)> callback_fn;
		typedef std::map <std::string, std::size_t> sample_name_map;
		
	protected:
		struct fsm
		{
			char	*p{nullptr};
			char	*pe{nullptr};
			char	*eof{nullptr};
		};
		
		template <typename> struct caller;
		template <typename> friend struct caller;
		
		template <int t_continue, int t_break>
		int check_max_field(vcf_field const field, int const target, callback_fn const &cb);

		void report_unexpected_character(char const *current_character, int const current_state);

	protected:
		fsm							m_fsm;
		transient_variant			m_current_variant;
		sample_name_map				m_sample_names;
		std::vector <char>			m_buffer;
		std::vector <format_field>	m_format;
		std::istream				*m_stream{nullptr};
		char						*m_line_start{nullptr};		// Current line start.
		char const					*m_start{0};				// Current string start.
		std::istream::pos_type		m_first_variant_offset{0};
		std::size_t					m_last_header_lineno{0};
		std::size_t					m_lineno{0};
		std::size_t					m_sample_idx{0};			// Current sample idx (1-based).
		std::size_t					m_idx{0};					// Current index in multi-part fields.
		std::size_t					m_len{0};
		std::size_t					m_pos{0};
		std::size_t					m_format_idx{0};
		std::size_t					m_integer{0};				// Currently read from the input.
		vcf_field					m_max_parsed_field{};
		bool						m_gt_is_phased{false};		// Is the current GT phased.
		bool						m_alt_is_complex{false};	// Is the current ALT “complex” (includes *).
	
	public:
		vcf_reader():
			m_buffer(128)
		{
		}
		
		vcf_reader(std::istream &stream):
			m_buffer(128),
			m_stream(&stream)
		{
		}
		
		void set_stream(std::istream &stream) { m_stream = &stream; }
		void read_header();
		void fill_buffer();
		void reset();
		bool parse(callback_fn &&callback);
		bool parse(callback_fn const &callback);
		
		std::size_t lineno() const { return m_lineno; }
		size_t sample_no(std::string const &sample_name) const;
		size_t sample_count() const { return m_sample_names.size(); }
		sample_name_map const &sample_names() const { return m_sample_names; }
		void set_parsed_fields(vcf_field max_field) { m_max_parsed_field = max_field; }
		
	protected:
		void skip_to_next_nl();
	};
}

#endif
