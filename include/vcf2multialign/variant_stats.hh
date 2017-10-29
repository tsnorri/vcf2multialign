/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_STATS_HH
#define VCF2MULTIALIGN_VARIANT_STATS_HH

#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/sequence_writer.hh>
#include <vcf2multialign/status_logger.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class error_logger;
	
	
	class variant_stats : public virtual sequence_writer_delegate
	{
	protected:
		variant_set							m_overlapping_alts;
		std::vector <skipped_sample>		m_skipped_samples;			// In current variant.
		std::map <uint8_t, sample_count>	m_counts_by_alt;			// In current variant.
		sample_count						m_non_ref_totals;			// In current variant.
		
	public:
		virtual class error_logger &error_logger() = 0;
		virtual class status_logger &status_logger() = 0;

		void handle_variant(variant_base const &var);
		
		// sequence_writer_delegate
		virtual void assigned_alt_to_sequence(std::size_t const alt_idx) override;
		virtual void found_overlapping_alt(
			std::size_t const lineno,
			uint8_t const alt_idx,
			std::size_t const sample_no,
			uint8_t const chr_idx
		) override;
		virtual void handled_alt(std::size_t const alt_idx) override;
		virtual void handled_haplotypes(variant_base const &var) override;
	};
}

#endif
