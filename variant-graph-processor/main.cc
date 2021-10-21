/*
 * Copyright (c) 2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/file_handling.hh>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;
namespace vgs	= vcf2multialign::variant_graphs;


namespace {
	
	constexpr inline std::uint64_t BRIDGE_NODE_POSITION_MASK{~std::uint64_t(0) >> 1};
	
	// Store the flag in the second argument to the highest bit.
	inline std::uint64_t make_bridge_node_record(std::uint64_t const pos, bool const has_only_ref_edges)
	{
		// Check that the highest bit is not in use in the position.
		libbio_always_assert_eq(pos, pos & BRIDGE_NODE_POSITION_MASK);
		
		std::uint64_t retval{has_only_ref_edges};
		retval <<= 63;
		retval |= pos;
		return retval;
	}


	void output_sample_names(vgs::variant_graph const &graph)
	{
		// Output the sample names.
		for (auto const &name : graph.sample_names())
			std::cout << name << '\n';
	}
	
	
	void output_bridge_nodes(vgs::variant_graph const &graph)
	{
		// Output the node positions in MSA co-ordinates.
		// We could just write the positions to std::cout but we would like to have the count, too.
		// The simplest way is to use a vector.
		
		std::vector <std::uint64_t> positions;
		vgs::variant_graph_walker walker(graph);
		
		walker.setup();
		std::size_t max_alt_edge_target{};
		while (vgs::variant_graph_walker::state::END != walker.advance())
		{
			auto const node(walker.node()); // 0-based.
			auto const &alt_edge_targets(walker.alt_edge_targets());
			
			// Check for parallel ALT edges.
			if (max_alt_edge_target <= node)
				positions.push_back(make_bridge_node_record(walker.ref_position(), ranges::empty(alt_edge_targets)));
			
			// Store the max. ALT edge target.
			for (auto const dst : alt_edge_targets)
				max_alt_edge_target = std::max(max_alt_edge_target, dst);
		}
		
		// Output.
		cereal::PortableBinaryOutputArchive archive(std::cout);
		archive(positions);
		std::cout << std::flush;
	}
}


int main(int argc, char **argv)
{
#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif

	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	vgs::variant_graph graph;
	
	{
		lb::file_istream input_graph_stream;
		
		// Open the graph file and read.
		lb::open_file_for_reading(args_info.variants_arg, input_graph_stream);
		cereal::PortableBinaryInputArchive iarchive(input_graph_stream);
		iarchive(graph);
	}
	
	if (args_info.bridge_nodes_given)
		output_bridge_nodes(graph);
	else if (args_info.sample_names_given)
		output_sample_names(graph);
	else
	{
		std::cerr << "No mode given.\n";
		std::exit(EXIT_FAILURE);
	}
	
	return EXIT_SUCCESS;
}
