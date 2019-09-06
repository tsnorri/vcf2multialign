/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_HAPLOTYPES_CHECK_OVERLAPPING_NON_NESTED_VARIANTS_HH
#define VCF2MULTIALIGN_HAPLOTYPES_CHECK_OVERLAPPING_NON_NESTED_VARIANTS_HH

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/haplotypes/error_logger.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	size_t check_overlapping_non_nested_variants(
		libbio::vcf_reader &reader,
		std::string const &chromosome_name,
		variant_set /* out */ &skipped_variants,
		error_logger &error_logger
	);
}

#endif
