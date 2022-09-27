/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <numeric>
#include "variant_filter.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void variant_filter::clear_queue()
	{
		// mnv_combiner needs the variants in this order.
		for (auto &var : m_msa_output_variants)
			m_next_handler->handle_variant_description(std::move(var));
		
		for (auto &var : m_vc_output_variants)
			m_next_handler->handle_variant_description(std::move(var));
		
		m_msa_output_variants.clear();
		m_vc_output_variants.clear();
	}
	
	
	void variant_filter::handle_variant_description(variant_description &&desc)
	{
		// Omit filtering for now, except for checking REF against ALT.
		if (desc.ref == desc.alt)
			desc.filters.emplace_back("ALT_EQ_TO_REF");
		
		if (0 == std::accumulate(desc.genotype.begin(), desc.genotype.end(), std::uint16_t(0)))
			desc.filters.emplace_back("GT_NOT_SET");
		
		libbio_assert_lte(m_current_pos, desc.position);
		if (m_current_pos < desc.position)
			clear_queue();
		
		m_current_pos = desc.position;
		switch (desc.origin)
		{
			case variant_origin::MSA:
				m_msa_output_variants.emplace_back(std::move(desc));
				break;
			
			case variant_origin::VC:
				m_vc_output_variants.emplace_back(std::move(desc));
				break;
		}
	}
	
	
	void variant_filter::finish()
	{
		clear_queue();
		m_next_handler->finish();
	}
}
