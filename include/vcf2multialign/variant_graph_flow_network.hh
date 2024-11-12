/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_FLOW_NETWORK_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_FLOW_NETWORK_HH

#include <climits>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vcf2multialign/variant_graph.hh>


namespace vcf2multialign::variant_graphs {

	// Adapts the variant graph by adding a source node and a sink node.
	// To make using the graph easier, REF edges are also added to the out-edge list.
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

		enum edge_property_mask : edge_type
		{
			special			= edge_type(1) << (CHAR_BIT * sizeof(edge_type) - 1),
			property_value	= special - 1
		};

		enum edge_property : edge_type
		{
			ref_edge					= edge_property_mask::special | 1,
			reverse_ref_edge			= edge_property_mask::special | 2,
			reverse_alt_edge			= edge_property_mask::special | 3,
			supplementary_edge			= edge_property_mask::special | 4,
			reverse_supplementary_edge	= edge_property_mask::special | 5
		};

		variant_graph const		&graph;
		node_vector				edge_sources;
		node_vector				edge_targets;
		edge_vector				edge_properties;			// ALT edge number in the original graph or the edge type.
		edge_vector				reverse_edges;				// Reverse edges
		edge_vector				out_edge_count_csum;		// Cumulative sum of out-edge counts by 1-based node number.

		explicit flow_network(variant_graph const &graph_):
			graph(graph_)
		{
		}

		void prepare();

		std::size_t node_count() const { return out_edge_count_csum.size() - 1; }
		std::size_t edge_count() const { return out_edge_count_csum.back(); }

		std::pair <edge_type, edge_type> out_edge_range(node_type const node) const { return {out_edge_count_csum[node], out_edge_count_csum[node + 1]}; }

		edge_type properties(edge_type edge) const { return edge_properties[edge]; }
	};
}

#endif
