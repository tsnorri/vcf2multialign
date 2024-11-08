/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/bits.hh>
#include <libbio/utility.hh>
#include <numeric>
#include <utility>
#include <vcf2multialign/transpose_matrix.hh>
#include <vcf2multialign/variant_graph.hh>
#include <vcf2multialign/variant_graph_flow_network.hh>
#include <vcf2multialign/variant_graph_flow_network_bgl_adapter.hh>

namespace lb	= libbio;


namespace {

	typedef vcf2multialign::variant_graphs::flow_network flow_network_type;


	struct edge_capacity_map
	{
		typedef flow_network_type::edge_type		key_type;
		typedef flow_network_type::capacity_type	value_type;

		flow_network_type const &flow_network;

		value_type operator[](key_type const key) const;
	};


	struct edge_weight_map
	{
		typedef flow_network_type::edge_type	key_type;
		typedef flow_network_type::weight_type	value_type;

		flow_network_type const &flow_network;

		value_type operator[](key_type const key) const;
	};


	auto edge_capacity_map::operator[](key_type const key) const -> value_type
	{
		if (0 == key)
			return flow_network.requested_flow;

		auto const is_reverse(1 == key % 2);
		if (is_reverse)
			return 0;

		auto const edge_idx(flow_network.variant_graph_edge_indices[key / 2]);
		auto const col(flow_network.graph.paths_by_edge_and_chrom_copy.column(edge_idx));
		auto const capacity(std::accumulate(col.word_begin(), col.word_end(), value_type(0), [](value_type acc, auto word) -> value_type {
			return acc + libbio::bits::count_bits_set(word);
		}));
		return capacity;
	}


	auto edge_weight_map::operator[](key_type const key) const -> value_type
	{
		auto const is_reverse(1 == key % 2);
		if (is_reverse)
			return 0;

		auto const edge_idx(key / 2);
		auto const edge_idx_(flow_network.variant_graph_edge_indices[edge_idx]);
		if (flow_network_type::EDGE_MAX == edge_idx_ || flow_network_type::REF_EDGE == edge_idx_)
			return 0;

		auto const src_idx(flow_network.edge_sources[edge_idx] - 1);
		auto const dst_idx(flow_network.edge_targets[edge_idx] - 1);

		auto const &ref_positions(flow_network.graph.reference_positions);
		value_type const ref_len(ref_positions[dst_idx] - ref_positions[src_idx]);
		value_type const alt_len(flow_network.graph.alt_edge_labels[edge_idx_].size());

		// Other options for applying weights to the edges include:
		// – Unit score for ALT edges, zero for REF edges
		// – Absolute value of ALT length minus REF length (works for indels, not for substitutions)
		// – Edit distance (difficult to calculate)
		return -1 * std::max(ref_len, alt_len);
	}


	inline edge_capacity_map::value_type get(
		edge_capacity_map const &map,
		edge_capacity_map::key_type const key
	)
	{
		return map[key];
	}


	inline edge_weight_map::value_type get(
		edge_weight_map const &map,
		edge_weight_map::key_type const key
	)
	{
		return map[key];
	}
}


namespace boost {

	template <>
	struct property_traits <edge_capacity_map>
	{
		typedef lvalue_property_map_tag				category;
		typedef edge_capacity_map::key_type			key_type;
		typedef edge_capacity_map::value_type		value_type;
		typedef edge_capacity_map::value_type		reference;
	};


	template <>
	struct property_traits <edge_weight_map>
	{
		typedef lvalue_property_map_tag				category;
		typedef edge_weight_map::key_type			key_type;
		typedef edge_weight_map::value_type			value_type;
		typedef edge_weight_map::value_type			reference;
	};
}


namespace vcf2multialign {

	// Apply a (very simple) algorithm to phase the variants in the given graph.
	// The algorithm currently supports one sample (not checked) and works as follows:
	// – The variant graph is first transformed into a flow network. A source and a sink node are added and the capacity of the single
	//   edge from the source node is set to the expected ploidy.
	// – The capacity of each REF edge is set to infinite and the weight to zero.
	// – The capacity of each ALT edge is set to the sum of the GT values that correspond to the edge and the weight is set to -||REF| - |ALT||.
	// – A minimum cost flow through the network is then calculated and edges are assigned to each chromosome copy based on positive flow.
	void phase(variant_graph &graph, std::uint16_t const ploidy)
	{
		variant_graphs::flow_network flow_network(graph, ploidy);
		edge_capacity_map const edge_capacities(flow_network);
		edge_weight_map const edge_weights(flow_network);
		boost::vector_property_map <variant_graphs::flow_network::node_type> vertex_predecessors(flow_network.node_count());
		boost::vector_property_map <variant_graphs::flow_network::capacity_type> edge_residual_capacities(flow_network.edge_count());
		boost::vector_property_map <boost::default_color_type> vertex_colors(flow_network.node_count());
		boost::vector_property_map <std::int64_t> vertex_distances(flow_network.node_count());

		// Return the flow for the given node.
		auto const flow_value([&edge_capacities, &edge_residual_capacities](variant_graphs::flow_network::edge_type const edge_idx){
			auto const capacity(edge_capacities[edge_idx]);
			auto const residual(edge_residual_capacities[edge_idx]);
			return capacity - residual;
		});

		// Decrease the flow through the given node by one unit.
		auto const decrease_flow([&edge_residual_capacities](variant_graphs::flow_network::edge_type const edge_idx){
			++edge_residual_capacities[edge_idx];
		});

		lb::log_time(std::cerr) << "Building the flow network…\n";
		flow_network.prepare();

		lb::log_time(std::cerr) << "Calculating the maximum flow…\n";
		boost::boykov_kolmogorov_max_flow(
			flow_network,
			0,
			flow_network.node_count() - 1,
			boost::capacity_map(edge_capacities)
			.residual_capacity_map(edge_residual_capacities)
			.predecessor_map(vertex_predecessors)
			.color_map(vertex_colors)
			.distance_map(vertex_distances)
		);

		lb::log_time(std::cerr) << "Calculating the minimum weight flow…\n";
		boost::cycle_canceling(
			flow_network,
			boost::residual_capacity_map(edge_residual_capacities)
			.weight_map(edge_weights)
			.predecessor_map(vertex_predecessors)
			.distance_inf(vertex_distances)
		);

		auto const calculated_flow(flow_value(0));
		if (calculated_flow < ploidy)
		{
			std::cerr << "ERROR: Unable to find " << ploidy << " paths through the variant graph; found " << calculated_flow << ".\n";
			std::exit(EXIT_FAILURE);
		}

		auto const &paths_by_edge_and_chrom_copy(graph.paths_by_edge_and_chrom_copy); // Edges on rows, chromosome copies (samples multiplied by ploidy) in columns.
		variant_graph::path_matrix new_paths_by_edge_and_chrom_copy(paths_by_edge_and_chrom_copy.number_of_rows(), ploidy);

		// We try to start from an arbitrary edge in order to distribute the variants more evenly to each chromosome copy.
		{
			std::uint64_t edge_idx{};
			for (std::uint16_t chr_idx{}; chr_idx < ploidy; ++chr_idx)
			{
				auto current_chr(new_paths_by_edge_and_chrom_copy.column(chr_idx));

				variant_graphs::flow_network::node_type node_idx{}; // Corresponds to the original graph.
				variant_graphs::flow_network::node_type const node_limit{graph.node_count()};
				while (node_idx < node_limit)
				{
					auto const out_edge_range(flow_network.out_edge_range(node_idx + 1)); // Take the source node into account.
					auto const out_edge_count(out_edge_range.second - out_edge_range.first);
					auto const current_edge_limit(edge_idx + out_edge_count);
					while (edge_idx < current_edge_limit)
					{
						auto const flow(flow_value(edge_idx % out_edge_count));
						if (flow)
						{
							decrease_flow(edge_idx);
							auto const edge_idx_(flow_network.variant_graph_edge_indices[edge_idx]);
							libbio_assert_neq(edge_idx_, variant_graphs::flow_network::EDGE_MAX);
							++edge_idx; // Not used before the next iteration, so we update here.

							if (variant_graphs::flow_network::REF_EDGE == edge_idx_)
								++node_idx;
							else
							{
								current_chr[edge_idx_] |= 1;
								node_idx = graph.alt_edge_targets[edge_idx_];
							}

							goto continue_traversal;
						}

						++edge_idx;
					}

					// Not found.
					libbio_fail("Unable to find a node with remaining flow.");

				continue_traversal:
					;
				}
			}
		}

		graph.paths_by_edge_and_chrom_copy = std::move(new_paths_by_edge_and_chrom_copy);
		graph.paths_by_chrom_copy_and_edge = transpose_matrix(graph.paths_by_edge_and_chrom_copy);
	}
}
