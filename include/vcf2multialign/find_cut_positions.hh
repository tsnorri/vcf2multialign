/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_FIND_CUT_POSITIONS_HH
#define VCF2MULTIALIGN_FIND_CUT_POSITIONS_HH

#include <limits>
#include <vcf2multialign/variant_graph.hh>
#include <vector>


namespace vcf2multialign {

	typedef variant_graph::ploidy_type cut_position_score_type;
	constexpr static inline auto const CUT_POSITION_SCORE_MAX{std::numeric_limits <cut_position_score_type>::max()};

	typedef std::vector <variant_graph::position_type> cut_position_vector;


	cut_position_score_type find_initial_cut_positions_lambda_min(
		variant_graph const &graph,
		variant_graph::edge_type const min_length,
		cut_position_vector &out_cut_positions,
		process_graph_delegate &delegate
	);
}

#endif
