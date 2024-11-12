/*
 * Copyright (c) 2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_FLOW_NETWORK_BGL_ADAPTER_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_FLOW_NETWORK_BGL_ADAPTER_HH

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <utility>
#include <vcf2multialign/variant_graph_flow_network.hh>

namespace vcf2multialign::variant_graphs::inline bgl {

	typedef vcf2multialign::variant_graphs::flow_network	flow_network;
}


namespace boost {

	template <>
	struct graph_traits <vcf2multialign::variant_graphs::bgl::flow_network>
	{
		struct traversal_category :
			public virtual vertex_list_graph_tag,
			public virtual edge_list_graph_tag,
			public virtual incidence_graph_tag
		{
		};

		typedef vcf2multialign::variant_graphs::bgl::flow_network	graph_type;
		typedef graph_type::node_type								vertex_descriptor;
		typedef graph_type::edge_type								edge_descriptor;

		typedef directed_tag										directed_category;
		typedef allow_parallel_edge_tag								edge_parallel_category;
		typedef graph_type::node_type								vertices_size_type;
		typedef graph_type::edge_type								edges_size_type;
		typedef edges_size_type										degree_size_type; // We could use some smaller type.

		typedef iterators::counting_iterator <
			graph_type::node_type
		>															vertex_iterator;

		typedef iterators::counting_iterator <
			graph_type::edge_type
		>															edge_iterator;

		typedef edge_iterator										out_edge_iterator;

		static vertex_descriptor null_vertex() { return graph_type::NODE_MAX; }
	};
}


namespace vcf2multialign::variant_graphs::inline bgl {

	struct vertex_index_map
	{
		typedef flow_network::node_type		key_type;
		typedef flow_network::node_type 	value_type;

		value_type operator[](key_type const key) const { return key; }
	};


	struct edge_reverse_map
	{
		typedef flow_network::edge_type		key_type;
		typedef flow_network::edge_type 	value_type;

		flow_network const &flow_network_;

		value_type operator[](key_type const edge) const { return flow_network_.reverse_edges[edge]; }
	};


	// VertexListGraph
	inline typename boost::graph_traits <flow_network>::vertices_size_type
	num_vertices(flow_network const &flow_network_)
	{
		return flow_network_.node_count();
	}


	inline std::pair <
		typename boost::graph_traits <flow_network>::vertex_iterator,
		typename boost::graph_traits <flow_network>::vertex_iterator
	>
	vertices(flow_network const &flow_network_)
	{
		return {{0}, {flow_network_.node_count()}};
	}


	// EdgeListGraph
	inline typename boost::graph_traits <flow_network>::edges_size_type
	num_edges(flow_network const &flow_network_)
	{
		return flow_network_.edge_count();
	}


	inline std::pair <
		typename boost::graph_traits <flow_network>::edge_iterator,
		typename boost::graph_traits <flow_network>::edge_iterator
	>
	edges(flow_network const &flow_network_)
	{
		return {{0}, {flow_network_.edge_count()}};
	}


	inline typename boost::graph_traits <flow_network>::vertex_descriptor
	source(
		typename boost::graph_traits <flow_network>::edge_descriptor const edge,
		flow_network const &flow_network_
	) {  return flow_network_.edge_sources[edge]; }


	inline typename boost::graph_traits <flow_network>::vertex_descriptor
	target(
		typename boost::graph_traits <flow_network>::edge_descriptor const edge,
		flow_network const &flow_network_
	) {  return flow_network_.edge_targets[edge]; }


	// IncidenceGraph
	inline std::pair <
		typename boost::graph_traits <flow_network>::out_edge_iterator,
		typename boost::graph_traits <flow_network>::out_edge_iterator
	>
	out_edges(
		typename boost::graph_traits <flow_network>::vertex_descriptor const node,
		flow_network const &flow_network_
	)
	{
		auto const range(flow_network_.out_edge_range(node));
		return {{range.first}, {range.second}};
	}


	inline typename boost::graph_traits <flow_network>::edges_size_type
	out_degree(
		typename boost::graph_traits <flow_network>::vertex_descriptor const node,
		flow_network const &flow_network_
	)
	{
		return (flow_network_.out_edge_count_csum[node + 1] - flow_network_.out_edge_count_csum[node]);
	}


	// Properties
	inline vertex_index_map get(boost::vertex_index_t, flow_network const &)
	{
		return {};
	}


	inline edge_reverse_map get(boost::edge_reverse_t, flow_network const &flow_network_)
	{
		return {flow_network_};
	}


	inline vertex_index_map::value_type get(vertex_index_map const &map, vertex_index_map::key_type const key)
	{
		return map[key];
	}


	inline edge_reverse_map::value_type get(edge_reverse_map const &map, edge_reverse_map::key_type const key)
	{
		return map[key];
	}
}


namespace boost {

	// Properties
	template <>
	struct property_map <vcf2multialign::variant_graphs::flow_network, vertex_index_t>
	{
		typedef vcf2multialign::variant_graphs::bgl::vertex_index_map const_type;
	};


	template <>
	struct property_map <vcf2multialign::variant_graphs::flow_network, edge_reverse_t>
	{
		typedef vcf2multialign::variant_graphs::bgl::edge_reverse_map const_type;
	};


	template <>
	struct property_traits <vcf2multialign::variant_graphs::bgl::vertex_index_map>
	{
		typedef lvalue_property_map_tag									category;
		typedef vcf2multialign::variant_graphs::bgl::vertex_index_map	map_type;
		typedef map_type::key_type										key_type;
		typedef map_type::value_type									value_type;
		typedef map_type::value_type									reference;
	};


	template <>
	struct property_traits <vcf2multialign::variant_graphs::bgl::edge_reverse_map>
	{
		typedef lvalue_property_map_tag									category;
		typedef vcf2multialign::variant_graphs::bgl::edge_reverse_map	map_type;
		typedef map_type::key_type										key_type;
		typedef map_type::value_type									value_type;
		typedef map_type::value_type									reference;
	};
}

#endif
