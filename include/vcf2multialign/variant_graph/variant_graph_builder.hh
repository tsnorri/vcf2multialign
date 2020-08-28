/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_VARIANT_GRAPH_BUILDER_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_VARIANT_GRAPH_BUILDER_HH

#include <string>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include <vector>


namespace vcf2multialign { namespace variant_graphs {
	
	struct builder_alt_edge
	{
		std::size_t	target_node{};
		std::string	label;
	
		builder_alt_edge() = default;
	
		builder_alt_edge(std::size_t const target_node_, std::string const &label_):
			target_node(target_node_),
			label(label_)
		{
		}
	};
	
	
	struct builder_node_description
	{
		struct subgraph_start_node_tag {};
		
		std::vector <builder_alt_edge>		alt_edges;
		std::vector <std::size_t>			paths_by_sample;
		std::vector <std::size_t>			edges_by_path;
		std::string							ref;
		std::size_t							pos{};
		bool								begins_subgraph{};
	
		builder_node_description() = default;
		
		// Create a node that begins a subgraph.
		builder_node_description(
			std::size_t pos_,
			std::string ref_,
			std::vector <builder_alt_edge> alt_edges_,
			std::vector <std::size_t> edges_by_path_,
			std::vector <std::size_t> paths_by_sample_,
			subgraph_start_node_tag const &
		):
			alt_edges(std::move(alt_edges_)),
			paths_by_sample(std::move(paths_by_sample_)),
			edges_by_path(std::move(edges_by_path_)),
			ref(std::move(ref_)),
			pos(pos_),
			begins_subgraph(true)
		{
		}
		
		// Create a node within a subgraph.
		builder_node_description(
			std::size_t pos_,
			std::string ref_,
			std::vector <builder_alt_edge> alt_edges_,
			std::vector <std::size_t> edges_by_path_
		):
			alt_edges(std::move(alt_edges_)),
			edges_by_path(std::move(edges_by_path_)),
			ref(std::move(ref_)),
			pos(pos_),
			begins_subgraph(false)
		{
		}
	};
	
	
	// Simple (to some extent) interface for building a variant graph (from nodes, not VCF records).
	class variant_graph_builder
	{
	public:
		variant_graph build_graph(
			std::vector <builder_node_description> const &nodes,
			std::vector <std::string> const &sample_names
		);
	};
}}

#endif
