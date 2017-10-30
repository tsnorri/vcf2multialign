/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_ALT_CHECKER_HH
#define VCF2MULTIALIGN_ALT_CHECKER_HH

#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant.hh>
#include <vector>


namespace vcf2multialign {

	class alt_checker
	{
	protected:
		std::vector <std::vector <uint8_t>>	m_valid_alts_by_lineno;
		variant_set const					*m_skipped_variants{nullptr};
		error_logger						*m_error_logger{};
		sv_handling							m_sv_handling_method{};
		std::size_t							m_max_alt_field_size{0};		// Need not be atomic b.c. set only when parsing the file.
		std::size_t							m_records_with_valid_alts{0};	// Need not be atomic b.c. set only when parsing the file.
		std::size_t							m_last_header_lineno{0};
		
	public:
		alt_checker() = default;
		
		alt_checker(
			std::size_t const lines,
			std::size_t const last_header_lineno,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants,
			error_logger &error_logger
		):
			m_valid_alts_by_lineno(lines),
			m_skipped_variants(&skipped_variants),
			m_error_logger(&error_logger),
			m_sv_handling_method(sv_handling_method),
			m_last_header_lineno(last_header_lineno)
		{
		}
		
		void check_variant(transient_variant const &var);
		void check_all_variants(vcf_reader &reader);
		
		std::size_t max_alt_field_size() const { return m_max_alt_field_size; }
		std::size_t	records_with_valid_alts() const { return m_records_with_valid_alts; }
		
		// variant_handler_delegate compatible, thread-safe.
		std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const
		{
			assert(m_last_header_lineno < lineno);
			return m_valid_alts_by_lineno.at(lineno - m_last_header_lineno - 1);
		}
		
		bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const
		{
			auto const &va(valid_alts(lineno));
			return (std::find(va.cbegin(), va.cend(), alt_idx) != va.cend());
		}
		
	protected:
		bool check_alt_seq(std::string_view const &alt) const;
	};
}

#endif
