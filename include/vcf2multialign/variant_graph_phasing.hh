/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_PHASING_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_PHASING_HH

#include <cstdint>
#include <vcf2multialign/variant_graph.hh>
#include <vcf2multialign/variant_graph_flow_network.hh>
#include <vector>


namespace vcf2multialign::variant_graphs::phasing {

	typedef vcf2multialign::variant_graphs::flow_network flow_network_type;


	struct edge_capacity_map
	{
		typedef flow_network_type::edge_type		key_type;
		typedef flow_network_type::capacity_type	value_type;

		flow_network_type const *flow_network{};
		value_type max_capacity{};

		value_type operator[](key_type const key) const;
		void output(flow_network_type const &) const;
	};


	struct edge_weight_map
	{
		typedef flow_network_type::edge_type	key_type;
		typedef flow_network_type::weight_type	value_type;

		flow_network_type const *flow_network{};

		value_type operator[](key_type const key) const;
		void output(flow_network_type const &) const;
	};
}


namespace vcf2multialign::variant_graphs {

	struct graph_phasing_delegate;


	class graph_phasing
	{
	public:
		typedef variant_graphs::flow_network		flow_network_type;
		typedef flow_network_type::edge_type		edge_type;
		typedef flow_network_type::capacity_type	edge_capacity_type;
		typedef flow_network_type::weight_type		edge_weight_type;

	private:
		phasing::edge_capacity_map					m_edge_capacities;
		phasing::edge_weight_map					m_edge_weights;
		std::vector <edge_capacity_type>			m_edge_residual_capacities;
		variant_graph								*m_graph{};

	public:
		explicit graph_phasing(variant_graph &graph):
			m_graph(&graph)
		{
		}

		bool phase(std::uint16_t const ploidy, graph_phasing_delegate &delegate);

		edge_capacity_type edge_capacity(edge_type edge) const { return m_edge_capacities[edge]; }
		edge_capacity_type edge_residual_capacity(edge_type edge) const { return m_edge_residual_capacities[edge]; }
		edge_capacity_type edge_flow(edge_type edge) const { return edge_capacity(edge) - edge_residual_capacity(edge); }
		edge_weight_type edge_weight(edge_type edge) const { return m_edge_weights[edge]; }

	private:
		void decrease_flow(edge_type const edge) { ++m_edge_residual_capacities[edge]; }
		void find_paths(flow_network_type const &flow_network, variant_graph::path_matrix &new_paths_by_edge_and_chrom_copy, std::uint16_t ploidy);
	};


	struct graph_phasing_delegate
	{
		typedef graph_phasing::flow_network_type	flow_network_type;

		virtual ~graph_phasing_delegate() {}
		virtual void graph_phasing_will_build_flow_network(graph_phasing const &) = 0;
		virtual void graph_phasing_will_calculate_maximum_flow(graph_phasing const &) = 0;
		virtual void graph_phasing_will_calculate_minimum_weight_flow(graph_phasing const &) = 0;
		virtual void graph_phasing_did_calculate_minimum_weigth_flow(graph_phasing const &, flow_network_type const &flow_network) = 0;
		virtual void graph_phasing_will_determine_paths(graph_phasing const &) = 0;

		virtual void graph_phasing_unable_to_match_ploidy(graph_phasing const &, std::uint16_t const ploidy, std::uint16_t const calculated_flow) = 0;
	};
}

#endif
