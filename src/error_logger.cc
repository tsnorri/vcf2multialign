/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/error_logger.hh>


namespace v2m = vcf2multialign;


namespace vcf2multialign {
	
	void error_logger::write_header()
	{
		m_output_stream
		<< "VARIANT_1_LINE\t"
		<< "VARIANT_2_LINE\t"
		<< "REF_POS\t"
		<< "SAMPLE\t"
		<< "CHR\t"
		<< "HANDLED_WITH_ALT\t"
		<< "TOTAL_WITH_ALT\t"
		<< "HANDLED_WITH_ANY_ALT\t"
		<< "TOTAL_WITH_ANY_ALT\t"
		<< "REASON\n";
	}

	
	void error_logger::log_conflicting_variants(std::size_t const line_1, std::size_t const line_2)
	{
		m_output_stream << line_1 << '\t' << line_2 << "\t\t\t\t\t\t\t\t" << "Conflicting variants\n";
	}
	
	
	void error_logger::log_ref_mismatch(std::size_t const lineno, std::size_t const diff_pos)
	{
		m_output_stream << lineno << "\t\t" << diff_pos
		<< "\t\t\t\t\t\t\tREF does not match the reference sequence (output anyway using the reference)\n";
	}

	
	void error_logger::log_overlapping_alternative(
		std::size_t const lineno,
		std::size_t const sample_no,
		std::size_t const chr_idx,
		sample_count const &alt_counts,
		sample_count const &non_ref_total_counts
	)
	{
		m_output_stream
		<< lineno << "\t\t\t"
		<< sample_no << '\t'
		<< chr_idx << '\t'
		<< alt_counts.handled_count << '\t'
		<< alt_counts.total_count << '\t'
		<< non_ref_total_counts.handled_count << '\t'
		<< non_ref_total_counts.total_count << '\t'
		<< "Overlapping alternative\n";
	}
}
