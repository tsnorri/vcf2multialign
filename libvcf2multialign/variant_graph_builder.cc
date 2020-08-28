/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/variant_graph/variant_graph_builder.hh>

namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	
	template <typename t_range>
	std::tuple <std::size_t, std::size_t> count_and_max(t_range &&range)
	{
		std::size_t count(0);
		std::size_t max(0);	// FIXME: could be the value type of the range.
		
		for (auto const &val : range)
		{
			++count;
			max = std::max(max, val);
		}
		return {count, max};
	}
}


namespace vcf2multialign { namespace variant_graphs {

	variant_graph variant_graph_builder::build_graph(
		std::vector <builder_node_description> const &nodes,
		std::vector <std::string> const &src_sample_names
	)
	{
		variant_graph retval;
		
		auto &alt_edge_labels(retval.alt_edge_labels());
		auto &alt_edge_targets(retval.alt_edge_targets());
		auto &sample_paths(retval.sample_paths());
		auto &path_edges(retval.path_edges());
		auto &sample_names(retval.sample_names());
		auto const &subgraph_start_positions(retval.subgraph_start_positions());
		std::size_t current_subgraph_idx(SIZE_MAX);
		std::size_t current_variant_idx_1(0); // Within the subgraph, 1-based.
		
		// Copy the sample names.
		sample_names = src_sample_names;
		
		// Handle the node descriptions.
		for (auto const &[i, node_desc] : rsv::enumerate(nodes))
		{
			auto const [node_idx, alt_edge_start] = retval.add_main_node(node_desc.pos, node_desc.alt_edges.size());
			libbio_assert_eq(node_idx, i);
			if (node_desc.begins_subgraph)
			{
				if (SIZE_MAX != current_subgraph_idx)
					retval.update_aligned_ref_positions_from_node(subgraph_start_positions[current_subgraph_idx]);
				
				// Count nodes with ALTs in the current subgraph.
				auto const [nodes_with_alts_tail, max_alts_tail] = count_and_max(
					nodes
					| rsv::drop_exactly(1 + i) // Remove the first 1 + i elements.
					| rsv::take_while([](auto const &nd){ return !nd.begins_subgraph; }) // Consider nodes until the next subgraph start.
					| rsv::filter([](auto const &nd){ return 0 < nd.alt_edges.size(); }) // Consider nodes that have ALT edges.
					| rsv::transform([](auto const &nd) -> std::size_t { return nd.alt_edges.size(); })
				);
				auto const nodes_with_alts((node_desc.alt_edges.size() ? 1 : 0) + nodes_with_alts_tail);
				auto const max_alts(std::max(max_alts_tail, node_desc.alt_edges.size()));
				
				// Add the subgraph.
				current_subgraph_idx = retval.add_subgraph(
					node_idx,
					node_desc.paths_by_sample.size(),
					node_desc.edges_by_path.size()
				);
				retval.setup_path_edges_for_current_subgraph(max_alts, nodes_with_alts, node_desc.edges_by_path.size());
				
				for (auto const [sample_idx, path_idx] : rsv::enumerate(node_desc.paths_by_sample))
					sample_paths[current_subgraph_idx][sample_idx] |= path_idx;
				
				current_variant_idx_1 = 0;
			}
			
			// Update variant count if the current node is a variant, i.e. has ALT edges.
			if (node_desc.alt_edges.size())
				++current_variant_idx_1;
			
			// Update the ALT edges
			for (auto const &[j, edge] : rsv::enumerate(node_desc.alt_edges))
			{
				alt_edge_targets[alt_edge_start + j] = edge.target_node;
				alt_edge_labels[alt_edge_start + j] = edge.label;
			}
			
			// Update the paths.
			if (SIZE_MAX != current_subgraph_idx)
			{
				libbio_assert_lt(0, current_variant_idx_1);
				for (auto const &[path_idx, edge_idx] : rsv::enumerate(node_desc.edges_by_path))
					path_edges[current_subgraph_idx](current_variant_idx_1 - 1, path_idx) |= edge_idx;
			}
		}
		
		// Update aligned positions.
		if (SIZE_MAX != current_subgraph_idx)
			retval.update_aligned_ref_positions_from_node(subgraph_start_positions[current_subgraph_idx]);
		
		return retval;
	}
}}
