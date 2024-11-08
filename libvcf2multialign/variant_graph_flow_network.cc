/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstddef>
#include <libbio/assert.hh>
#include <numeric>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>
#include <vcf2multialign/variant_graph.hh>
#include <vcf2multialign/variant_graph_flow_network.hh>

namespace rsv	= ranges::views;


namespace vcf2multialign::variant_graphs {

	void flow_network::prepare()
	{
		edge_targets.clear();
		edge_sources.clear();
		reverse_edges.clear();
		edge_properties.clear();
		out_edge_count_csum.clear();

		edge_targets.resize(2 * (1 + graph.node_count() + graph.edge_count()), NODE_MAX); // Edges for source and sink nodes, too.
		edge_sources.resize(edge_targets.size(), NODE_MAX);
		reverse_edges.resize(edge_targets.size(), EDGE_MAX);
		edge_properties.resize(edge_targets.size(), edge_property_mask::property_value);
		out_edge_count_csum.resize(graph.node_count() + 3, 0); // 1-based with source and sink nodes

		edge_type edge_idx{}; // Zero for the first edge from source to the first node.

		auto const add_edge_with_reverse([this, &edge_idx](node_type src, node_type dst, edge_type edge_properties_, edge_type reverse_edge_properties){
			edge_sources[edge_idx] = src;
			edge_targets[edge_idx] = dst;
			reverse_edges[edge_idx] = edge_idx + 1;
			edge_properties[edge_idx] = edge_properties_;
			++edge_idx;

			// Reverse
			edge_sources[edge_idx] = dst;
			edge_targets[edge_idx] = src;
			reverse_edges[edge_idx] = edge_idx - 1;
			edge_properties[edge_idx] = reverse_edge_properties;
			++edge_idx;

			++out_edge_count_csum[1 + src];
			++out_edge_count_csum[1 + dst];
		});

		// Source node
		add_edge_with_reverse(0, 1, flow_network::edge_property::supplementary_edge, flow_network::edge_property::reverse_supplementary_edge);

		{
			variant_graph_walker walker(graph);
			while (walker.advance())
			{
				add_edge_with_reverse(walker.node() + 1, walker.node() + 2, flow_network::edge_property::ref_edge, flow_network::edge_property::reverse_ref_edge);

				for (auto const [idx, target_node] : rsv::enumerate(walker.alt_edge_targets()))
					add_edge_with_reverse(walker.node() + 1, target_node + 1, walker.alt_edge_base() + idx, flow_network::edge_property::reverse_alt_edge);
			}
		}

		// Sink node.
		edge_properties[edge_idx - 2] = flow_network::edge_property::supplementary_edge;
		edge_properties[edge_idx - 1] = flow_network::edge_property::reverse_supplementary_edge;

		libbio_assert_eq(edge_targets.size(), edge_idx);

		// Sort the edges by the source node.
		{
			edge_vector permutation(edge_sources.size());
			edge_vector inverse(edge_sources.size());
			std::iota(permutation.begin(), permutation.end(), 0);
			ranges::sort(rsv::zip(edge_sources, edge_targets, reverse_edges, edge_properties, permutation));
			for (auto const [ii, jj] : rsv::enumerate(permutation))
				inverse[jj] = ii;

			// Apply the inverse permutation to the reverse edges.
			for (auto &re : reverse_edges)
				re = inverse[re];
		}

		// Calculate the cumulative sum.
		{
			auto const size(out_edge_count_csum.size());
			for (std::size_t i(1); i < size; ++i)
				out_edge_count_csum[i] += out_edge_count_csum[i - 1];
		}
	}
}
