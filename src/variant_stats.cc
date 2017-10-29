/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/variant_stats.hh>


namespace vcf2multialign {
	
	void variant_stats::assigned_alt_to_sequence(std::size_t const alt_idx)
	{
		++m_non_ref_totals.handled_count;
		++m_counts_by_alt[alt_idx].handled_count;
	}
	
	
	void variant_stats::found_overlapping_alt(
		std::size_t const lineno,
		uint8_t const alt_idx,
		std::size_t const sample_no,
		uint8_t const chr_idx
	)
	{
		if (m_overlapping_alts.insert(lineno).second)
		{
			status_logger().log([lineno, sample_no, chr_idx](){
				std::cerr << "Overlapping alternatives on line " << lineno
					<< " for sample " << sample_no << ':' << (int) chr_idx
					<< " (and possibly others); skipping when needed." << std::endl;
			});
		}
	
		if (error_logger().is_logging_errors())
			m_skipped_samples.emplace_back(sample_no, alt_idx, chr_idx);
	}
	
	
	void variant_stats::handle_variant(variant_base const &var)
	{
		m_skipped_samples.clear();
		m_counts_by_alt.clear();
		m_non_ref_totals.reset();
	}
	
	
	void variant_stats::handled_alt(std::size_t const alt_idx)
	{
		++m_non_ref_totals.total_count;
		++m_counts_by_alt[alt_idx].total_count;
	}
	
	
	void variant_stats::handled_haplotypes(variant_base const &var)
	{
		// Report errors if needed.
		auto &el(error_logger());
		if (el.is_logging_errors())
		{
			auto const lineno(var.lineno());
			for (auto const &s : m_skipped_samples)
				el.log_overlapping_alternative(lineno, s.sample_no, s.chr_idx, m_counts_by_alt[s.alt_idx], m_non_ref_totals);
		}
	}
}
