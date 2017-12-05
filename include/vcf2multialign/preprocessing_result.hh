/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESSING_RESULT_HH
#define VCF2MULTIALIGN_PREPROCESSING_RESULT_HH

#include <vcf2multialign/alt_checker.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {

	struct preprocessing_result
	{
		class alt_checker	alt_checker;
		vector_type			reference;
		ploidy_map			ploidy;
		variant_set			skipped_variants;
		subgraph_map		subgraph_starting_points;	// Variant line numbers by character pointers.
	};
}

#endif
