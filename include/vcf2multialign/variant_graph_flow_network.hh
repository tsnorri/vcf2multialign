/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_FLOW_NETWORK_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_FLOW_NETWORK_HH

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vcf2multialign/variant_graph.hh>


namespace vcf2multialign::variant_graphs {

	struct flow_network
	{
		typedef variant_graph::node_type	node_type;
		typedef variant_graph::edge_type	edge_type;
		typedef variant_graph::node_vector	node_vector;
		typedef variant_graph::edge_vector	edge_vector;
		typedef std::int32_t				capacity_type;
		typedef std::int32_t				weight_type;

		constexpr static inline node_type const NODE_MAX{variant_graph::NODE_MAX};
		constexpr static inline edge_type const EDGE_MAX{variant_graph::EDGE_MAX};
		constexpr static inline edge_type const REF_EDGE{EDGE_MAX - 1};

		variant_graph const		&graph;
		node_vector				edge_sources;
		node_vector				edge_targets;
		edge_vector				variant_graph_edge_indices;	// From the flow network edges
		edge_vector				flow_network_edge_indices;	// From the original graph edges
		edge_vector				edge_count_csum;			// Cumulative sum of ALT edge counts by 1-based node number.

		capacity_type			requested_flow{};

		flow_network(variant_graph const &graph_, capacity_type const requested_flow_):
			graph(graph_),
			requested_flow(requested_flow_)
		{
		}

		void prepare();

		std::size_t node_count() const { return edge_count_csum.size() - 1; }
		std::size_t edge_count() const { return edge_count_csum.back(); }

		std::pair <edge_type, edge_type> out_edge_range(node_type const node) const { return {edge_count_csum[node], edge_count_csum[node + 1]}; }
	};
}

#endif
