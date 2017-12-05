/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CHECK_OVERLAPPING_NON_NESTED_VARIANTS_HH
#define VCF2MULTIALIGN_CHECK_OVERLAPPING_NON_NESTED_VARIANTS_HH

#include <vcf2multialign/logger.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {
	std::size_t check_overlapping_non_nested_variants(
		vcf_reader &reader,
		sv_handling const sv_handling_method,
		variant_set /* out */ &skipped_variants,
		struct logger &logger
	);
}

#endif
