/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/bits.hh>
#include <vcf2multialign/preprocess/variant_graph.hh>


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
		m_ref_positions.reserve(1 + expected_count);
		m_aligned_ref_positions.reserve(1 + expected_count);
		m_alt_edge_count_csum.reserve(1 + expected_count);
		
		m_alt_edge_targets.reserve(2 * expected_count);
		m_alt_edge_labels.reserve(2 * expected_count);
	}
	
	
	std::size_t variant_graph::add_subgraph(std::size_t const node_idx, std::size_t const sample_count, std::size_t const variant_count, std::size_t const path_count)
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
		auto &path_edge_matrix(m_path_edges.emplace_back(path_count, variant_count, path_bits));
		
		for (auto &word : sample_path_vec.word_range())
			word = 0;
		for (auto &word : path_edge_matrix.word_range())
			word = 0;
		
		return subgraph_idx;
	}
	
	
	std::tuple <std::size_t, std::size_t, bool> variant_graph::add_main_node(std::size_t const ref_pos, std::size_t const alt_edge_count)
	{
		auto const prev_ref_pos(m_ref_positions.back());
		
		libbio_always_assert_lte(prev_ref_pos, ref_pos);
		libbio_always_assert(!m_alt_edge_count_csum.empty());
		
		auto const prev_alt_edge_count(m_alt_edge_count_csum.back());
		auto const new_alt_edge_count(prev_alt_edge_count + alt_edge_count);
		
		// Check if the requested node already exists.
		bool did_create(false);
		if (prev_ref_pos == ref_pos)
		{
			// Adjust the ALT edge count.
			m_alt_edge_count_csum.back() = new_alt_edge_count;
		}
		else
		{
			// Otherwise create.
			did_create = true;
			m_ref_positions.emplace_back(ref_pos);
			m_aligned_ref_positions.resize(m_ref_positions.size(), 0);
			m_alt_edge_count_csum.emplace_back(new_alt_edge_count);
		}
		
		m_alt_edge_targets.resize(new_alt_edge_count, 0);
		m_alt_edge_labels.resize(new_alt_edge_count);
		
		// m_ref_positions’s indexing is 1-based but we would like to return the 0-based index.
		return std::tuple <std::size_t, std::size_t, bool>(m_ref_positions.size() - 2, prev_alt_edge_count, did_create);
	}
	
	
	void variant_graph::connect_alt_edges(std::size_t const src_idx, std::size_t const dst_idx)
	{
		auto const start(m_alt_edge_count_csum[src_idx]);
		auto const end(m_alt_edge_count_csum[1 + src_idx]);
		for (std::size_t i(start); i < end; ++i)
			m_alt_edge_targets[i] = dst_idx;
	}
}
