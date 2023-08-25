/*
 * Copyright (c) 2019–2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_MSA_COMBINER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_MSA_COMBINER_HH

#include <array>
#include <vcf2multialign/types.hh>
#include <vector>
#include "mnv_combiner.hh"
#include "msa_data_source.hh"
#include "overlap_counter.hh"
#include "variant_deduplicator.hh"
#include "variant_filter.hh"
#include "variant_writer.hh"
#include "vcf_record_generator.hh"
#include "types.hh"


namespace vcf2multialign {
	
	class msa_combiner : public segmentation_handler
	{
		friend msa_data_source;

	public:
		struct combining_statistics
		{
			std::size_t	combined_variants{};
			std::size_t	alt_eq_to_ref{};
		};
		
	protected:
		msa_data_source						m_data_source;
		
		// Output handling
		variant_writer						m_variant_writer;		// output_handler
		mnv_combiner						m_mnv_combiner;			// output_handler
		variant_filter						m_variant_filter;		// output_handler
		variant_deduplicator				m_variant_deduplicator;
		
		std::size_t							m_alt_matches_ref{};
		std::uint16_t						m_ploidy{1};
	
	public:
		msa_combiner() = default;
		
		msa_combiner(
			std::string output_chr_id,
			std::uint16_t const ploidy,
			std::ostream &os,
			bool const should_output_msa_variants,
			bool const should_log_status
		):
			m_data_source(ploidy, should_log_status),
			m_variant_writer(os, std::move(output_chr_id), should_output_msa_variants),
			m_mnv_combiner(m_variant_writer, ploidy),
			m_variant_filter(m_mnv_combiner),
			m_variant_deduplicator(m_variant_filter),
			m_ploidy(ploidy)
		{
		}
		
		// Entry point.
		combining_statistics process_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen);
		
	protected:
		void process_one_segment_msa(
			aligned_segment const &seg,
			std::int32_t overlap_count,
			overlap_counter::const_iterator overlap_it,
			overlap_counter::const_iterator const overlap_end
		);
		
		std::size_t process_variants_msa(
			std::size_t const max_alt_pos,
			aligned_segment_vector::const_iterator seg_it,
			aligned_segment_vector::const_iterator const seg_end,
			class overlap_counter const &overlap_counter
		);
		
		void process_variant_in_range(
			variant_record const &var,
			aligned_segment_vector::const_iterator aligned_segment_begin,
			aligned_segment_vector::const_iterator const aligned_segment_end,
			std::int32_t const max_overlaps
		);

		void process_variants(msa_segmentation &segmentation) override;
	};
}

#endif