/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/error_logger.hh>


namespace v2m = vcf2multialign;


namespace vcf2multialign {
	
	void error_logger::write_header()
	{
		m_output_stream << "VARIANT_1_LINE\tVARIANT_2_LINE\tSAMPLE\tCHR\tHANDLED_NON_REF_SAMPLE_COUNT\tREASON\n";
	}

	
	void error_logger::log_conflicting_variants(std::size_t const line_1, std::size_t const line_2)
	{
		m_output_stream << line_1 << '\t' << line_2 << "\t\t\t\t" << "Conflicting variants\n";
	}

	
	void error_logger::log_overlapping_alternative(
		std::size_t const lineno,
		std::size_t const sample_no,
		uint8_t const chr_idx,
		std::size_t const handled_non_ref_samples
	)
	{
		m_output_stream << lineno << "\t\t" << sample_no << '\t' << chr_idx << '\t' << handled_non_ref_samples << "\tOverlapping alternative\n";
	}
}
