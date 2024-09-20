/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/state.hh>

namespace ml	= libbio::memory_logger;


namespace vcf2multialign {
	
	char const *to_chars(state const state_)
	{
		switch (state_)
		{
			case state::default_state: 						return "default_state";
			case state::build_variant_graph: 				return "build_variant_graph";
			case state::output_haplotypes: 					return "output_haplotypes";
			case state::output_founder_sequences_greedy:	return "output_founder_sequences_greedy";
			case state::find_cut_positions: 				return "find_cut_positions";
			case state::find_matchings: 					return "find_matchings";
			case state::state_limit:						return "state_limit";
		}
		
		return "unknown";
	}
}
