/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_FIND_FIRST_MATCHING_VARIANT_HH
#define VCF2MULTIALIGN_UTILITY_FIND_FIRST_MATCHING_VARIANT_HH

#include <libbio/vcf/vcf_reader.hh>


namespace vcf2multialign {
	bool find_first_matching_variant(libbio::vcf::reader &reader, std::string_view const chr_id);
	bool find_variant_by_line_number(libbio::vcf::reader &reader, std::size_t const lineno);
}

#endif
