/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_COMBINE_MSA_HH
#define VCF2MULTIALIGN_COMBINE_MSA_COMBINE_MSA_HH

#include "msa_combiner.hh"
#include "vcf_record_generator.hh"

namespace vcf2multialign {
	
	void combine_msa(
		char const *ref_path,
		char const *alt_path,
		char const *variants_path,
		char const *output_chr,
		std::uint16_t const ploidy,
		std::ostream &os
	);
}

#endif
