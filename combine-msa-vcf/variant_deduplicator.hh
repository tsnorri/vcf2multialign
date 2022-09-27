/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VARIANT_DEDUPLICATOR_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VARIANT_DEDUPLICATOR_HH

#include <vector>
#include "variant_description.hh"


namespace vcf2multialign {
	
	struct output_handler;
	
	
	class variant_deduplicator
	{
	protected:
		std::vector <variant_description>	m_output_variants;
		output_handler						*m_next_handler{};
		
	public:
		variant_deduplicator() = default;
		variant_deduplicator(output_handler &next_handler):
			m_next_handler(&next_handler)
		{
		}
		
		variant_description &handle_variant_description(variant_description &&desc) { return m_output_variants.emplace_back(std::move(desc)); }

		std::size_t size() const { return m_output_variants.size(); }
		
		void merge_output_variants(std::size_t const partition_point);
		void filter_processed_variants_and_output(std::size_t const min_unhandled_ref_pos);
		void finish();
	};
}

#endif
