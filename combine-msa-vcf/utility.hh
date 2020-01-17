/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_UTILITY_HH
#define VCF2MULTIALIGN_COMBINE_MSA_UTILITY_HH

#include <libbio/vcf/variant.hh>


namespace vcf2multialign {
	
	std::pair <std::uint16_t, std::uint16_t> count_set_genotype_values(libbio::variant const &var, std::uint16_t const alt_idx);
}

#endif
