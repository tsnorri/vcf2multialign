/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_FIND_OPTIMAL_CUT_POSITIONS_HH
#define VCF2MULTIALIGN_PREPROCESS_FIND_OPTIMAL_CUT_POSITIONS_HH

#include <string>
#include <vector>


namespace vcf2multialign {
	
	void find_optimal_cut_positions(
		char const *reference_path,
		char const *variants_path,
		char const *reference_seq_name,
		char const *chr_name,
		std::vector <std::string> const &field_names_for_filter_if_set,
		std::size_t const minimum_subgraph_distance,
		bool const should_overwrite_files
	);
}

#endif
