/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstddef>
#include <libbio/assert.hh>
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
		variant_graph_edge_indices.clear();
		flow_network_edge_indices.clear();
		edge_count_csum.clear();

		edge_targets.resize(2 + graph.node_count() + graph.edge_count(), NODE_MAX); // Edges for source and sink nodes, too.
		edge_sources.resize(edge_targets.size(), NODE_MAX);
		variant_graph_edge_indices.resize(edge_targets.size(), EDGE_MAX);
		flow_network_edge_indices.resize(graph.edge_count());
		edge_count_csum.resize(graph.node_count() + 3, 0); // 1-based with source and sink nodes

		edge_type edge_idx{}; // Zero for the first edge from source to the first node.

		auto const add_edge_with_reverse([this, &edge_idx](node_type src, node_type dst, edge_type original_edge){
			edge_sources[edge_idx] = src;
			edge_targets[edge_idx] = dst;
			variant_graph_edge_indices[edge_idx] = original_edge;
			flow_network_edge_indices[original_edge] = edge_idx;
			++edge_idx;

			++edge_count_csum[1 + src];
			++edge_count_csum[1 + dst];
		});

		// Source node
		add_edge_with_reverse(0, 1, EDGE_MAX);

		{
			variant_graph_walker walker(graph);
			while (walker.advance())
			{
				add_edge_with_reverse(walker.node() + 1, walker.node() + 2, REF_EDGE);

				for (auto const [idx, target_node] : rsv::enumerate(walker.alt_edge_targets()))
					add_edge_with_reverse(walker.node() + 1, target_node + 1, walker.alt_edge_base() + idx);
			}
		}

		// Sink node
		add_edge_with_reverse(graph.node_count(), graph.node_count() + 1, EDGE_MAX);

		libbio_assert_eq(edge_targets.size(), edge_idx);

		ranges::sort(rsv::zip(edge_sources, edge_targets));

		// Calculate the cumulative sum.
		{
			auto const size(edge_count_csum.size());
			for (std::size_t i(1); i < size; ++i)
				edge_count_csum[i] += edge_count_csum[i - 1];
		}
	}
}
