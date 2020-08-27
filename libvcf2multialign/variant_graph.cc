/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/bits.hh>
#include <vcf2multialign/variant_graph/variant_graph.hh>


namespace lb	= libbio;


namespace vcf2multialign {
	
	void variant_graph::clear()
	{
		m_ref_positions.clear();
		m_aligned_ref_positions.clear();
		m_alt_edge_targets.clear();
		m_alt_edge_count_csum.clear();
		m_subgraph_start_positions.clear();
		m_alt_edge_labels.clear();
		m_sample_names.clear();
		m_sample_paths.clear();
		m_path_edges.clear();
		
		// Add entries for node zero and its ALT edge count.
		m_alt_edge_count_csum.resize(2, 0);
		
		// Both of these are 1-based, i.e. node zero’s position is stored at index 1.
		m_ref_positions.resize(2, 0);
		m_aligned_ref_positions.resize(2, 0);
	}
	
	
	void variant_graph::reserve_memory_for_nodes(std::size_t const expected_count)
	{
		// FIXME: call me.
		m_ref_positions.reserve(1 + expected_count);
		m_aligned_ref_positions.reserve(1 + expected_count);
		m_alt_edge_count_csum.reserve(1 + expected_count);
		
		m_alt_edge_targets.reserve(2 * expected_count);
		m_alt_edge_labels.reserve(2 * expected_count);
	}
	
	
	std::size_t variant_graph::add_subgraph(std::size_t const node_idx, std::size_t const sample_count, std::size_t const path_count)
	{
		libbio_assert_lt(0, path_count);
		libbio_assert(0 == m_subgraph_start_positions.size() || m_subgraph_start_positions.back() < node_idx);
		m_subgraph_start_positions.emplace_back(node_idx);
		auto const subgraph_idx(m_subgraph_start_positions.size() - 1);
		
		// Check that the counts are in sync.
		libbio_assert_eq(m_sample_paths.size(), subgraph_idx);
		libbio_assert_eq(m_path_edges.size(), subgraph_idx);
		auto const path_bits(lb::bits::highest_bit_set(path_count - 1) ?: 1);
		auto &sample_path_vec(m_sample_paths.emplace_back(sample_count, path_bits));
		
		for (auto &word : sample_path_vec.word_range())
			word = 0;
		
		m_max_paths_in_subgraph = std::max(m_max_paths_in_subgraph, path_count);
		
		return subgraph_idx;
	}
	
	
	void variant_graph::setup_path_edges_for_current_subgraph(std::size_t const max_node_alt_count, std::size_t const variant_count, std::size_t const path_count)
	{
		auto const edge_bits(lb::bits::highest_bit_set(max_node_alt_count) ?: 1);
		auto &path_edge_matrix(m_path_edges.emplace_back(variant_count, path_count, edge_bits));
		for (auto &word : path_edge_matrix.word_range())
			word = 0;
	}
	
	
	std::tuple <std::size_t, std::size_t> variant_graph::add_main_node(std::size_t const ref_pos, std::size_t const alt_edge_count)
	{
		auto const prev_ref_pos(m_ref_positions.back());
		
		libbio_always_assert_lt(prev_ref_pos, ref_pos);
		libbio_always_assert(!m_alt_edge_count_csum.empty());
		
		auto const prev_alt_edge_count(m_alt_edge_count_csum.back());
		auto const new_alt_edge_count(prev_alt_edge_count + alt_edge_count);
		
		m_ref_positions.emplace_back(ref_pos);
		m_aligned_ref_positions.resize(m_ref_positions.size(), 0);
		m_alt_edge_count_csum.emplace_back(new_alt_edge_count);
		
		m_alt_edge_targets.resize(new_alt_edge_count, 0);
		m_alt_edge_labels.resize(new_alt_edge_count);
		
		return {m_ref_positions.size() - 2, prev_alt_edge_count};
	}


	std::tuple <std::size_t, std::size_t> variant_graph::add_main_node_if_needed(std::size_t const ref_pos, std::size_t const alt_edge_count)
	{
		auto const prev_ref_pos(m_ref_positions.back());
		if (prev_ref_pos != ref_pos)
			return add_main_node(ref_pos, alt_edge_count);
		else
		{
			auto const ref_position_count(m_ref_positions.size());
			libbio_assert_lte(2, ref_position_count);
			auto const node_idx(ref_position_count - 2);
			
			auto const prev_alt_edge_count(m_alt_edge_count_csum.back());
			auto const new_alt_edge_count(prev_alt_edge_count + alt_edge_count);
			
			m_alt_edge_count_csum.back() = new_alt_edge_count;
			m_alt_edge_targets.resize(new_alt_edge_count, 0);
			m_alt_edge_labels.resize(new_alt_edge_count);
			
			return {m_ref_positions.size() - 2, prev_alt_edge_count};
		}
	}
	
	
	void variant_graph::update_aligned_ref_positions_from_node(std::size_t const first_node_idx)
	{
		// Calculate the aligned positions.
		// Handle nodes pairwise and at the same time update the aligned positions of the edge targets.
		// This works b.c. calculating the aligned position of the destination node of an (ALT) edge
		// based on that edge only requires the aligned position of the source node and the weight of
		// the edge to be known.
		
		libbio_assert_eq(m_ref_positions.size(), m_aligned_ref_positions.size());
		auto const total_node_count(m_ref_positions.size() - 1); // m_ref_positions uses 1-based indexing.
		libbio_assert_lt(1 + first_node_idx, total_node_count); // There should be at least two nodes (position and end of one variant).
		auto const node_pair_count(total_node_count - first_node_idx); // Start s.t. lhs is one node outside the current processed range.
		for (std::size_t i(0); i < node_pair_count; ++i)
		{
			// Determine the values for the current handled pair of indices.
			// Position indices are 1-based.
			auto const lhs_node_idx(first_node_idx + i - 1);
			auto const rhs_node_idx(1 + lhs_node_idx);
			auto const lhs_ref(m_ref_positions[1 + lhs_node_idx]);
			auto const lhs_aln(m_aligned_ref_positions[1 + lhs_node_idx]);
			auto const rhs_ref(m_ref_positions[1 + rhs_node_idx]);
			auto &rhs_aln(m_aligned_ref_positions[1 + rhs_node_idx]);
		
			// Update the aligned positions of the ALT edge targets.
			auto alt_idx(m_alt_edge_count_csum[lhs_node_idx]);
			auto const alt_limit(m_alt_edge_count_csum[rhs_node_idx]);
			while (alt_idx < alt_limit)
			{
				auto const target_node_idx(m_alt_edge_targets[alt_idx]);
				auto const label(m_alt_edge_labels[alt_idx]);
				auto &target_aln_pos(m_aligned_ref_positions[1 + target_node_idx]);
				target_aln_pos = std::max(target_aln_pos, lhs_aln + label.size());
				++alt_idx;
			}
			
			// Update the rhs aligned position.
			auto const ref_len(rhs_ref - lhs_ref);
			rhs_aln = std::max(rhs_aln, lhs_aln + ref_len);
		}
	}
	
	
	auto variant_graph_walker::advance() -> state
	{
		++m_node_1;
		if (m_graph->m_ref_positions.size() == m_node_1)
			return state::END;
		
		return state::NODE;
	}
	
	
	auto variant_graph_walker::advance_and_track_subgraph() -> state
	{
		++m_node_1;
		if (m_node_1 == m_next_subgraph_start_1)
		{
			++m_subgraph;
			
			auto const &subgraph_start_positions(m_graph->m_subgraph_start_positions);
			m_next_subgraph_start_1 = (
				m_subgraph < subgraph_start_positions.size()
				? 1 + subgraph_start_positions[m_subgraph]
				: SIZE_MAX
			);
			return state::SUBGRAPH_START_NODE;
		}
		
		if (m_graph->m_ref_positions.size() == m_node_1)
			return state::END;
		
		return state::NODE;
	}
}
