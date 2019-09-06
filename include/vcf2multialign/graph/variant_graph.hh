/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GRAPH_VARIANT_GRAPH_HH
#define VCF2MULTIALIGN_GRAPH_VARIANT_GRAPH_HH

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <libbio/int_vector.hh>
#include <libbio/int_vector/cereal_serialization.hh>
#include <libbio/int_matrix.hh>
#include <libbio/int_matrix/cereal_serialization.hh>
#include <tuple>
#include <vector>


namespace vcf2multialign {
	
	class variant_graph
	{
	public:
		typedef std::vector <std::size_t>	position_vector;
		typedef std::vector <std::string>	string_vector;
		typedef libbio::int_vector <0>		sample_path_vector;
		typedef libbio::int_matrix <0>		path_edge_matrix;
		
	protected:
		position_vector							m_ref_positions;			// REF positions (0-based) by node number. We use 1-based indexing in order to make summing easier.
		position_vector							m_aligned_ref_positions;	// Aligned REF positions by node number. We use 1-based indexing in order to make summing easier.
		position_vector							m_alt_edge_targets;			// ALT edge target nodes by edge number as a concatenated vector.
		position_vector							m_alt_edge_count_csum;		// Cumulative sum of ALT edges by node number.
		position_vector							m_subgraph_start_positions;	// Subgraph starting node numbers by subgraph number.
		string_vector							m_alt_edge_labels;			// ALT edge labels by edge number.
		string_vector							m_sample_names;				// Sample names by sample number.
		std::vector <sample_path_vector>		m_sample_paths;				// Sample path numbers by sample and subgraph number. (Subgraph count is not known at the start, so use a vector for them.)
		std::vector <path_edge_matrix>			m_path_edges;				// Edge numbers (0 for REF edge, 1 for first ALT edge etc.) by path, variant and subgraph number.
		std::size_t								m_max_paths_in_subgraph{};
		
		// Nodes that do not represent variant positions also have zero ALT edges.
		
	public:
		position_vector &ref_positions() { return m_ref_positions; }
		position_vector &aligned_ref_positions() { return m_aligned_ref_positions; }
		position_vector &alt_edge_targets() { return m_alt_edge_targets; }
		position_vector &alt_edge_count_csum() { return m_alt_edge_count_csum; }
		position_vector &subgraph_start_positions() { return m_subgraph_start_positions; }
		string_vector &alt_edge_labels() { return m_alt_edge_labels; }
		string_vector &sample_names() { return m_sample_names; }
		std::vector <libbio::int_vector <0>> &sample_paths() { return m_sample_paths; }
		std::vector <libbio::int_matrix <0>> &path_edges() { return m_path_edges; }
		
		position_vector const &ref_positions() const { return m_ref_positions; }
		position_vector const &aligned_ref_positions() const { return m_aligned_ref_positions; }
		position_vector const &alt_edge_targets() const { return m_alt_edge_targets; }
		position_vector const &alt_edge_count_csum() const { return m_alt_edge_count_csum; }
		position_vector const &subgraph_start_positions() const { return m_subgraph_start_positions; }
		string_vector const &alt_edge_labels() const { return m_alt_edge_labels; }
		string_vector const &sample_names() const { return m_sample_names; }
		std::vector <libbio::int_vector <0>> const &sample_paths() const { return m_sample_paths; }
		std::vector <libbio::int_matrix <0>> const &path_edges() const { return m_path_edges; }
		std::size_t max_paths_in_subgraph() const { return m_max_paths_in_subgraph; }
		
		void clear();
		void reserve_memory_for_nodes(std::size_t const expected_count);
		std::size_t add_subgraph(std::size_t const node_idx, std::size_t const sample_count, std::size_t const variant_count, std::size_t const path_count);
		std::tuple <std::size_t, std::size_t, bool> add_main_node(std::size_t const ref_pos, std::size_t const alt_edge_count);
		void connect_alt_edges(std::size_t const src_idx, std::size_t const dst_idx);
		
		// For Cereal.
		template <typename t_archive>
		void serialize(t_archive &archive, std::uint32_t const version);
	};
	
	
	template <typename t_archive>
	void variant_graph::serialize(t_archive &archive, std::uint32_t const version)
	{
		// Ignore the version for now.
		archive(
			m_ref_positions,
			m_aligned_ref_positions,
			m_alt_edge_targets,
			m_alt_edge_count_csum,
			m_subgraph_start_positions,
			m_alt_edge_labels,
			m_sample_names,
			m_sample_paths,
			m_path_edges,
			m_max_paths_in_subgraph
		);
	}
}

#endif
