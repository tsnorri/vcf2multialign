/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_FIND_SUBGRAPH_STARTING_POINTS_HH
#define VCF2MULTIALIGN_FIND_SUBGRAPH_STARTING_POINTS_HH

#include <vcf2multialign/types.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {
	void find_subgraph_starting_points(
		vcf_reader &reader,
		variant_set const &skipped_variants,
		variant_set /* out */ &subgraph_starting_points
	);
}

#endif
