/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_ERROR_LOGGER_HH
#define VCF2MULTIALIGN_ERROR_LOGGER_HH

#include <libbio/file_handling.hh>
#include <libbio/variant.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class error_logger
	{
	protected:
		libbio::file_ostream	m_output_stream;
		
	public:
		libbio::file_ostream &output_stream() { return m_output_stream; }
		void flush() { if (is_logging_errors()) m_output_stream.flush(); }
		bool is_logging_errors() const { return m_output_stream.is_open(); }
		void write_header();
		
		void log_no_supported_alts(std::size_t const line);
		
		void log_skipped_structural_variant(std::size_t const line, std::size_t const alt_idx, libbio::sv_type const svt);
		
		void log_invalid_alt_seq(std::size_t const line, std::size_t const alt_idx, std::string const &alt);
		
		void log_conflicting_variants(std::size_t const line_1, std::size_t const line_2);
		
		void log_ref_mismatch(std::size_t const lineno, std::size_t const diff_pos);
		
		void log_overlapping_alternative(
			std::size_t const lineno,
			std::size_t const sample_no,
			std::size_t const chr_idx,
			sample_count const &alt_counts,
			sample_count const &non_ref_total_counts
		);
	};
}

#endif
