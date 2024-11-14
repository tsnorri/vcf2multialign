/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/bits.hh>
#include <libbio/utility.hh>
#include <numeric>
#include <ostream>
#include <utility>
#include <vcf2multialign/transpose_matrix.hh>
#include <vcf2multialign/variant_graph.hh>
#include <vcf2multialign/variant_graph_flow_network.hh>
#include <vcf2multialign/variant_graph_flow_network_bgl_adapter.hh>
#include <vcf2multialign/variant_graph_phasing.hh>


namespace vcf2multialign::variant_graphs::phasing {

	auto edge_capacity_map::operator[](key_type const edge_idx) const -> value_type
	{
		auto const properties(flow_network->edge_properties[edge_idx]);
		switch (properties)
		{
			case flow_network_type::edge_property::ref_edge:
			case flow_network_type::edge_property::supplementary_edge:
				return max_capacity;

			case flow_network_type::edge_property::reverse_ref_edge:
			case flow_network_type::edge_property::reverse_supplementary_edge:
			case flow_network_type::edge_property::reverse_alt_edge:
				return 0;

			default: // ALT edge index
			{
				libbio_assert(! (flow_network_type::edge_property_mask::special & properties));
				auto const col(flow_network->graph.paths_by_edge_and_chrom_copy.column(properties));
				auto const capacity(std::accumulate(col.word_begin(), col.word_end(), value_type(0), [](value_type acc, auto word) -> value_type {
					return acc + libbio::bits::count_bits_set(word);
				}));
				return capacity;
			}
		}
	}


	auto edge_weight_map::operator[](key_type const edge_idx) const -> value_type
	{
		auto const properties(flow_network->edge_properties[edge_idx]);
		switch (properties)
		{
			case flow_network_type::edge_property::ref_edge:
			case flow_network_type::edge_property::reverse_ref_edge:
			case flow_network_type::edge_property::supplementary_edge:
			case flow_network_type::edge_property::reverse_supplementary_edge:
				return 0;

			case flow_network_type::edge_property::reverse_alt_edge:
				return -1 * operator[](flow_network->reverse_edges[edge_idx]);

			default: // ALT edge index
			{
				libbio_assert(! (flow_network_type::edge_property_mask::special & properties));

				auto const src_idx(flow_network->edge_sources[edge_idx] - 1);
				auto const dst_idx(flow_network->edge_targets[edge_idx] - 1);

				auto const &ref_positions(flow_network->graph.reference_positions);
				value_type const ref_len(ref_positions[dst_idx] - ref_positions[src_idx]);
				value_type const alt_len(flow_network->graph.alt_edge_labels[properties].size());

				// Other options for applying weights to the edges include:
				// – Unit score for ALT edges, zero for REF edges
				// – Absolute value of ALT length minus REF length (works for indels, not for substitutions)
				// – Edit distance (difficult to calculate)
				return -1 * std::max(ref_len, alt_len);
			}
		}
	}


	void edge_capacity_map::output(flow_network_type const &flow_network) const
	{
		std::cerr << "Edge capacities:\n";
		auto const edge_limit(flow_network.edge_count());
		for (flow_network_type::edge_type edge{}; edge < edge_limit; ++edge)
			std::cerr << "[" << edge << "]:\t" << this->operator[](edge) << '\n';
	}


	void edge_weight_map::output(flow_network_type const &flow_network) const
	{
		std::cerr << "Edge weights:\n";
		auto const edge_limit(flow_network.edge_count());
		for (flow_network_type::edge_type edge{}; edge < edge_limit; ++edge)
			std::cerr << "[" << edge << "]:\t" << this->operator[](edge) << '\n';
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


namespace {

	typedef vcf2multialign::variant_graphs::flow_network				flow_network_type;
	typedef vcf2multialign::variant_graphs::phasing::edge_capacity_map	edge_capacity_map;
	typedef vcf2multialign::variant_graphs::phasing::edge_weight_map	edge_weight_map;


	// Fix an issue in boost::cycle_canceling.
	// Using flow_network_type const & as the type of the first parameter seems to cause more problems.
	template <
		typename t_edge_residual_capacities,
		typename t_edge_weights,
		typename t_vertex_predecessors,
		typename t_vertex_distances
	>
	void cycle_canceling_(
		flow_network_type const &flow_network,
		t_edge_residual_capacities &edge_residual_capacities,
		t_edge_weights const &edge_weights,
		t_vertex_predecessors const &vertex_predecessors,
		t_vertex_distances const &vertex_distances

	)
	{
		auto const &reverse_edges([&]{
			using namespace boost;
			auto const &retval(get(edge_reverse, flow_network));
			return retval;
		}());

		boost::cycle_canceling(
			flow_network,
			edge_weights,
			reverse_edges,
			edge_residual_capacities,
			vertex_predecessors,
			vertex_distances
		);
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


namespace vcf2multialign::variant_graphs {

	void graph_phasing::find_paths(flow_network_type const &flow_network, variant_graph::path_matrix &new_paths_by_edge_and_chrom_copy, std::uint16_t const ploidy)
	{
		// We try to start from an arbitrary edge in order to distribute the variants more evenly to each chromosome copy.
		auto const &graph(flow_network.graph);
		std::uint64_t edge_idx{};
		libbio_assert(m_graph->node_count());
		for (std::uint16_t chr_idx{}; chr_idx < ploidy; ++chr_idx)
		{
			auto current_chr(new_paths_by_edge_and_chrom_copy.column(chr_idx));

			flow_network_type::node_type node_idx{}; // Corresponds to the original graph.
			flow_network_type::node_type const node_limit{m_graph->node_count() - 1}; // Do not process the out-edges of the last node (since they point to the sink).
			while (node_idx < node_limit)
			{
				auto const out_edge_range(flow_network.out_edge_range(node_idx + 1)); // Take the source node into account.
				auto const out_edge_count(out_edge_range.second - out_edge_range.first);
				auto const current_edge_limit(edge_idx + out_edge_count);
				while (edge_idx < current_edge_limit)
				{
					auto const edge_idx_(edge_idx % out_edge_count + out_edge_range.first);
					++edge_idx; // Not used before the next iteration, so we update here.

					auto const flow(edge_flow(edge_idx_));
					if (0 < flow)
					{
						auto const properties(flow_network.edge_properties[edge_idx_]);
						switch (properties)
						{
							case flow_network_type::edge_property::reverse_ref_edge:
							case flow_network_type::edge_property::reverse_alt_edge:
							case flow_network_type::edge_property::supplementary_edge:
							case flow_network_type::edge_property::reverse_supplementary_edge:
								libbio_fail("There should be no need to handle any reverse or supplementary forward edges while finding paths");
							case flow_network_type::edge_property::ref_edge:
								++node_idx;
								break;
							default: // ALT edge
							{
								libbio_assert_lt(1 + node_idx, graph.alt_edge_count_csum.size());
								auto const alt_edge_base(graph.alt_edge_count_csum[node_idx]);
								libbio_assert_lte(alt_edge_base, properties);
								libbio_assert_lt(properties, graph.alt_edge_count_csum[1 + node_idx]);
								current_chr[properties - alt_edge_base] |= 1;
								node_idx = graph.alt_edge_targets[properties];
								break;
							}
						}

						decrease_flow(edge_idx_);
						goto continue_traversal;
					}
				}

				// Not found.
				libbio_fail("Unable to find a node with remaining flow.");

			continue_traversal:
				;
			}
		}
	}


	// Apply a (very simple) algorithm to phase the variants in the given graph.
	// The algorithm currently supports one sample (not checked) and works as follows:
	// – The variant graph is first transformed into a flow network. A source and a sink node are added and the capacity of the single
	//   edge from the source node is set to the expected ploidy.
	// – The capacity of each REF edge is set to infinite and the weight to zero.
	// – The capacity of each ALT edge is set to the sum of the GT values that correspond to the edge and the weight is set to -||REF| - |ALT||.
	// – A minimum cost flow through the network is then calculated and edges are assigned to each chromosome copy based on positive flow.
	bool graph_phasing::phase(std::uint16_t const ploidy, graph_phasing_delegate &delegate)
	{
		delegate.graph_phasing_will_build_flow_network(*this);
		variant_graphs::flow_network flow_network(*m_graph);
		flow_network.prepare();

		m_edge_capacities = edge_capacity_map(&flow_network, ploidy);
		m_edge_weights = edge_weight_map(&flow_network);
		m_edge_residual_capacities.clear();
		m_edge_residual_capacities.resize(flow_network.edge_count(), 0);

		boost::vector_property_map <variant_graphs::flow_network::node_type> vertex_predecessors(flow_network.node_count());
		boost::vector_property_map <boost::default_color_type> vertex_colors(flow_network.node_count());
		boost::vector_property_map <std::int64_t> vertex_distances(flow_network.node_count());
		auto const &vertex_indices([&]{
			auto const &retval(get(boost::vertex_index, flow_network));
			return retval;
		}());
		auto edge_residual_capacities(boost::make_iterator_property_map(m_edge_residual_capacities.begin(), vertex_indices));

		delegate.graph_phasing_will_calculate_maximum_flow(*this);
		auto const calculated_flow(boost::boykov_kolmogorov_max_flow(
			flow_network,
			0,
			flow_network.node_count() - 1,
			boost::capacity_map(m_edge_capacities)
			.residual_capacity_map(edge_residual_capacities)
			.color_map(vertex_colors)
			.predecessor_map(vertex_predecessors)
			.distance_map(vertex_distances)
		));

		if (calculated_flow != ploidy)
		{
			delegate.graph_phasing_unable_to_match_ploidy(*this, ploidy, calculated_flow);
			return false;
		}

		delegate.graph_phasing_will_calculate_minimum_weight_flow(*this);
		cycle_canceling_(
			flow_network,
			edge_residual_capacities,
			m_edge_weights,
			vertex_predecessors,
			vertex_distances
		);
		delegate.graph_phasing_did_calculate_minimum_weigth_flow(*this, flow_network);

		{
			delegate.graph_phasing_will_determine_paths(*this);
			constexpr std::size_t const path_matrix_row_col_divisor{64}; // Make sure we can transpose the matrix with the 8×8 operation.
			std::size_t const path_matrix_cols(path_matrix_row_col_divisor * std::ceil(1.0 * ploidy / path_matrix_row_col_divisor));
			auto const &paths_by_edge_and_chrom_copy(m_graph->paths_by_edge_and_chrom_copy); // Edges on rows, chromosome copies (samples multiplied by ploidy) in columns.
			variant_graph::path_matrix new_paths_by_edge_and_chrom_copy(paths_by_edge_and_chrom_copy.number_of_rows(), path_matrix_cols);
			find_paths(flow_network, new_paths_by_edge_and_chrom_copy, ploidy);
			m_graph->paths_by_edge_and_chrom_copy = std::move(new_paths_by_edge_and_chrom_copy);
			m_graph->paths_by_chrom_copy_and_edge = transpose_matrix(m_graph->paths_by_edge_and_chrom_copy);
		}

		return true;
	}
}
