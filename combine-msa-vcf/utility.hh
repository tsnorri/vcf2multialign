/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_UTILITY_HH
#define VCF2MULTIALIGN_COMBINE_MSA_UTILITY_HH

#include <libbio/assert.hh>
#include <libbio/vcf/variant.hh>
#include <range/v3/all.hpp>
#include "algorithms.hh"
#include "types.hh"


namespace vcf2multialign {
	
	std::pair <std::uint16_t, std::uint16_t> count_set_genotype_values(libbio::variant const &var, std::uint16_t const alt_idx);
	void split_mixed_segment(aligned_segment &src, aligned_segment_vector &dst);
	
	
	inline void check_segment(aligned_segment const &seg)
	{
		namespace rsv	= ranges::view;
		
		// Check that for MISMATCH, every (aligned) character does not match.
		libbio_assert(
			segment_type::MISMATCH != seg.type ||
			ranges::all_of(
				rsv::zip(seg.ref.string, seg.alt.string),
				[](auto const &tup){
					return std::get <0>(tup) != std::get <1>(tup);
				}
			)
		);
	}
}

#endif
