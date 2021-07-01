/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_VARIANT_GRAPH_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_VARIANT_GRAPH_HH

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <libbio/int_vector.hh>
#include <libbio/int_vector/cereal_serialization.hh>
#include <libbio/int_matrix.hh>
#include <libbio/int_matrix/cereal_serialization.hh>
#include <tuple>
#include <vcf2multialign/types.hh>
#include <vector>


namespace vcf2multialign { namespace variant_graphs {
	
	class variant_graph
	{
		friend class variant_graph_walker;
		
	public:
		struct alt_edge
		{
			std::string const	&label;
			std::size_t			target_node{};
			
			alt_edge(std::size_t const target_node_, std::string const &label_):
				label(label_),
				target_node(target_node_)
			{
			}
		};
		
	public:
		typedef std::vector <std::size_t>		position_vector;
		typedef std::vector <std::string>		string_vector;
		typedef libbio::int_vector <0>			sample_path_vector;
		typedef libbio::int_matrix <0>			path_edge_matrix;
		
	protected:
		position_vector							m_ref_positions;			// REF positions (0-based) by node number. We use 1-based indexing in order to make summing easier.
		position_vector							m_aligned_ref_positions;	// Aligned REF positions by node number. We use 1-based indexing in order to make summing easier.
		position_vector							m_alt_edge_targets;			// ALT edge target nodes by edge number as a concatenated vector.
		position_vector							m_alt_edge_count_csum;		// Cumulative sum of ALT edges by node number.
		position_vector							m_subgraph_start_positions;	// Subgraph starting node numbers by subgraph number.
		string_vector							m_alt_edge_labels;			// ALT edge labels by edge number.
		string_vector							m_sample_names;				// Sample names by sample number.
		std::vector <sample_path_vector>		m_sample_paths;				// Sample path numbers by sample and subgraph number. (Subgraph count is not known at the start, so use a vector for them.)
		std::vector <path_edge_matrix>			m_path_edges;				// Edge numbers (0 for REF edge, 1 for first ALT edge etc.) by path (column), variant (row) and subgraph number.
		std::size_t								m_max_paths_in_subgraph{};
		
		// Nodes that do not represent variant positions also have zero ALT edges.
		// Variants correspond to nodes that do have ALT edges.
		
	public:
		std::size_t node_count() const { return m_ref_positions.size() - 1; }
		
		inline std::size_t ref_position_for_node(std::size_t const node_idx) const;
		inline std::size_t aligned_position_for_node(std::size_t const node_idx) const;
		inline std::string_view ref_label(std::size_t const lhs_node, std::size_t const rhs_node, vector_type const &ref) const;
		inline std::string_view ref_label(std::size_t const lhs_node, vector_type const &ref) const { return ref_label(lhs_node, 1 + lhs_node, ref); }
		inline auto alt_edge_labels(std::size_t const node_idx) const;
		inline auto alt_edge_targets(std::size_t const node_idx) const;
		inline auto alt_edges(std::size_t const node_idx) const;
		
		position_vector &ref_positions() { return m_ref_positions; }
		position_vector &aligned_ref_positions() { return m_aligned_ref_positions; }
		position_vector &alt_edge_targets() { return m_alt_edge_targets; }
		position_vector &alt_edge_count_csum() { return m_alt_edge_count_csum; }
		position_vector &subgraph_start_positions() { return m_subgraph_start_positions; }
		string_vector &alt_edge_labels() { return m_alt_edge_labels; }
		string_vector &sample_names() { return m_sample_names; }
		std::vector <sample_path_vector> &sample_paths() { return m_sample_paths; }
		std::vector <path_edge_matrix> &path_edges() { return m_path_edges; }
		
		position_vector const &ref_positions() const { return m_ref_positions; }
		position_vector const &aligned_ref_positions() const { return m_aligned_ref_positions; }
		position_vector const &alt_edge_targets() const { return m_alt_edge_targets; }
		position_vector const &alt_edge_count_csum() const { return m_alt_edge_count_csum; }
		position_vector const &subgraph_start_positions() const { return m_subgraph_start_positions; }
		string_vector const &alt_edge_labels() const { return m_alt_edge_labels; }
		string_vector const &sample_names() const { return m_sample_names; }
		std::vector <sample_path_vector> const &sample_paths() const { return m_sample_paths; }
		std::vector <path_edge_matrix> const &path_edges() const { return m_path_edges; }
		std::size_t subgraph_count() const { return m_subgraph_start_positions.size(); }
		std::size_t max_paths_in_subgraph() const { return m_max_paths_in_subgraph; }
		
		void clear();
		void reserve_memory_for_nodes(std::size_t const expected_count);
		std::tuple <std::size_t, std::size_t> add_main_node(std::size_t const ref_pos, std::size_t const alt_edge_count);
		std::tuple <std::size_t, std::size_t> add_main_node_if_needed(std::size_t const ref_pos, std::size_t const alt_edge_count);
		std::size_t add_subgraph(std::size_t const node_idx, std::size_t const sample_count, std::size_t const path_count);
		void setup_path_edges_for_current_subgraph(std::size_t const max_node_alt_count, std::size_t const variant_count, std::size_t const path_count);
		void update_aligned_ref_positions_from_node(std::size_t const first_node_idx);
		
		inline bool operator==(variant_graph const &other) const;
		
		// For Cereal.
		template <typename t_archive>
		void serialize(t_archive &archive, std::uint32_t const version);
	};
	
	
	bool variant_graph::operator==(variant_graph const &other) const
	{
		return
			m_ref_positions				== other.m_ref_positions &&
			m_aligned_ref_positions		== other.m_aligned_ref_positions &&
			m_alt_edge_targets			== other.m_alt_edge_targets &&
			m_alt_edge_count_csum		== other.m_alt_edge_count_csum &&
			m_subgraph_start_positions	== other.m_subgraph_start_positions &&
			m_alt_edge_labels			== other.m_alt_edge_labels &&
			m_sample_names				== other.m_sample_names &&
			m_sample_paths				== other.m_sample_paths &&
			m_path_edges				== other.m_path_edges &&
			m_max_paths_in_subgraph		== other.m_max_paths_in_subgraph;
	}
	
	
	auto variant_graph::alt_edge_labels(std::size_t const node_idx) const
	{
		libbio_assert_lt(1 + node_idx, m_alt_edge_count_csum.size());
		auto const start_idx(m_alt_edge_count_csum[node_idx]);
		auto const end_idx(m_alt_edge_count_csum[1 + node_idx]);
		return m_alt_edge_labels | ranges::view::slice(start_idx, end_idx);
	}
	
	
	auto variant_graph::alt_edge_targets(std::size_t const node_idx) const
	{
		libbio_assert_lt(1 + node_idx, m_alt_edge_count_csum.size());
		auto const start_idx(m_alt_edge_count_csum[node_idx]);
		auto const end_idx(m_alt_edge_count_csum[1 + node_idx]);
		return m_alt_edge_targets | ranges::view::slice(start_idx, end_idx);
	}
	
	
	auto variant_graph::alt_edges(std::size_t const node_idx) const
	{
		libbio_assert_lt(1 + node_idx, m_alt_edge_count_csum.size());
		auto const start_idx(m_alt_edge_count_csum[node_idx]);
		auto const end_idx(m_alt_edge_count_csum[1 + node_idx]);
		return
			ranges::view::zip(
				m_alt_edge_targets | ranges::view::slice(start_idx, end_idx),
				m_alt_edge_labels | ranges::view::slice(start_idx, end_idx)
			)
			| ranges::view::transform([](auto const &tup){
				auto const &[target_node, label] = tup;
				return alt_edge(target_node, label);
			});
	}
	
	
	class variant_graph_walker
	{
	public:
		enum state
		{
			NODE,
			SUBGRAPH_START_NODE,
			END
		};
		
	protected:
		vector_type const	*m_reference{};
		variant_graph const	*m_graph{};
		std::size_t			m_node_1{};							// 1-based.
		std::size_t			m_subgraph{};
		std::size_t			m_next_subgraph_start_1{SIZE_MAX};	// 1-based.
		
	public:
		variant_graph_walker() = default;

		variant_graph_walker(variant_graph const &graph):
			m_graph(&graph)
		{
		}
		
		variant_graph_walker(variant_graph const &graph, vector_type const &reference):
			m_reference(&reference),
			m_graph(&graph)
		{
		}
		
		void setup() { m_next_subgraph_start_1 = 1 + m_graph->m_subgraph_start_positions.front(); }
		state advance();
		state advance_and_track_subgraph();
		variant_graph const &graph() const { return *m_graph; }
		std::size_t node() const { return m_node_1 - 1; }
		std::size_t subgraph() const { return m_subgraph; }
		std::size_t ref_position() const { return m_graph->ref_position_for_node(m_node_1 - 1); }
		std::size_t aligned_position() const { return m_graph->aligned_position_for_node(m_node_1 - 1); }
		std::string_view ref_label() const { return ref_label_(m_node_1); }
		std::string_view ref_label(std::size_t const rhs_node) const { libbio_assert_lte(m_node_1 - 1, rhs_node); return ref_label_(rhs_node); }
		std::size_t ref_length() const { return ref_length_(m_node_1); }
		std::size_t ref_length(std::size_t const rhs_node) const { libbio_assert_lte(m_node_1 - 1, rhs_node); return ref_length_(rhs_node); }
		std::size_t aligned_length() const { return aligned_length_(m_node_1); }
		std::size_t aligned_length(std::size_t const rhs_node) const { libbio_assert_lte(m_node_1 - 1, rhs_node); return aligned_length_(rhs_node); }
		auto alt_edge_labels() const { return m_graph->alt_edge_labels(m_node_1 - 1); }
		auto alt_edge_targets() const { return m_graph->alt_edge_targets(m_node_1 - 1); }
		auto alt_edges() const { return m_graph->alt_edges(m_node_1 - 1); }
		
	protected:
		std::string_view ref_label_(std::size_t const rhs_node) const { return m_graph->ref_label(m_node_1 - 1, rhs_node, *m_reference); }
		inline std::size_t ref_length_(std::size_t const rhs_node) const;
		inline std::size_t aligned_length_(std::size_t const rhs_node) const;
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
	
	
	std::size_t variant_graph::ref_position_for_node(std::size_t const node_idx) const
	{
		libbio_assert_lt(1 + node_idx, m_ref_positions.size());
		return m_ref_positions[1 + node_idx];
	}
	
	
	std::size_t variant_graph::aligned_position_for_node(std::size_t const node_idx) const
	{
		libbio_assert_lt(1 + node_idx, m_aligned_ref_positions.size());
		return m_aligned_ref_positions[1 + node_idx];
	}
	
	
	std::string_view variant_graph::ref_label(std::size_t const lhs_node, std::size_t const rhs_node, vector_type const &ref) const
	{
		libbio_assert_lt(1 + lhs_node, m_ref_positions.size());
		libbio_assert_lt(1 + rhs_node, m_ref_positions.size());
		auto const lhs_ref_pos(m_ref_positions[1 + lhs_node]);
		auto const rhs_ref_pos(m_ref_positions[1 + rhs_node]);
		return std::string_view(ref.data() + lhs_ref_pos, rhs_ref_pos - lhs_ref_pos);
	}


	std::size_t variant_graph_walker::ref_length_(std::size_t const rhs_node) const
	{
		if (m_graph->m_ref_positions.size() <= 1 + rhs_node)
			return 0; // FIXME: add one element to m_ref_positions so that the check is not needed.
		return m_graph->ref_position_for_node(rhs_node) - m_graph->ref_position_for_node(m_node_1 - 1);
	}


	std::size_t variant_graph_walker::aligned_length_(std::size_t const rhs_node) const
	{
		if (m_graph->m_aligned_ref_positions.size() <= 1 + rhs_node)
			return 0; // FIXME: add one element to m_ref_positions so that the check is not needed.
		return m_graph->aligned_position_for_node(rhs_node) - m_graph->aligned_position_for_node(m_node_1 - 1);
	}
}}

#endif
