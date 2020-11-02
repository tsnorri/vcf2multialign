/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <numeric>
#include <range/v3/all.hpp>
#include "output_handler.hh"
#include "variant_filter.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void variant_filter::merge_output_variants(std::size_t const partition_point)
	{
		// Merge the partitions of sorted variants.
		std::inplace_merge(
			m_output_variants.begin(),
			m_output_variants.begin() + partition_point,
			m_output_variants.end(),
			[](auto const &lhs, auto const &rhs){
				return lhs.position < rhs.position;
			}
		);
		libbio_assert(
			std::is_sorted(
				m_output_variants.begin(),
				m_output_variants.end(),
				[](auto const &lhs, auto const &rhs){ return lhs.position < rhs.position; }
			)
		);
	}
	
	
	void variant_filter::filter_processed_variants_and_output(std::size_t const min_unhandled_ref_pos)
	{
		// Omit filtering for now, except for checking REF against ALT.
		libbio_assert(
			std::is_sorted(m_output_variants.begin(), m_output_variants.end(), [](auto const &lhs, auto const &rhs){
				return lhs.position < rhs.position;
			})
		);
		std::vector <std::string> filters;
		auto var_it(m_output_variants.begin());
		auto const var_end(m_output_variants.end());
		while (var_it != var_end)
		{
			auto &desc(*var_it);
			if (min_unhandled_ref_pos <= desc.position)
				break;
			
			if (desc.ref == desc.alt)
				desc.filters.emplace_back("ALT_EQ_TO_REF");
			
			if (0 == std::accumulate(desc.genotype.begin(), desc.genotype.end(), std::uint16_t(0)))
				desc.filters.emplace_back("GT_NOT_SET");
			
			m_next_handler->handle_variant_description(std::move(*var_it));
			++var_it;
		}
		
		m_output_variants.erase(m_output_variants.begin(), var_it);
	}
}
