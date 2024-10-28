/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_STATE_HH
#define VCF2MULTIALIGN_STATE_HH

#include <cstdint>
#include <libbio/log_memory_usage.hh>


namespace vcf2multialign {

	// For the memory logger
	enum class state : std::uint64_t
	{
		default_state = 0,
		build_variant_graph,
		output_haplotypes,
		output_founder_sequences_greedy,
		find_cut_positions,
		find_matchings,
		state_limit
	};


	char const *to_chars(state const state_);

	typedef libbio::memory_logger::header_writer_delegate_ <state> memory_logger_header_writer_delegate;
}

#endif
