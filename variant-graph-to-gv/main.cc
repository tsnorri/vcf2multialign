/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <range/v3/view/sliding.hpp>
#include <unistd.h>
#include <vcf2multialign/preprocess/variant_graph.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	void output_for_gv(
		char const *reference_path,
		char const *input_graph_path,
		char const *output_graph_path,
		char const *reference_seq_name,
		bool const should_overwrite_files
	)
	{
		v2m::vector_type reference;
		v2m::variant_graph graph;
		lb::file_ostream output_graph_stream;

		{
			lb::file_istream input_graph_stream;
			
			// Open the files.
			{
				auto const mode(lb::make_writing_open_mode({
					lb::writing_open_mode::CREATE,
					(should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
				}));
				lb::open_file_for_writing(output_graph_path, output_graph_stream, mode);
			}
			
			lb::open_file_for_reading(input_graph_path, input_graph_stream);
			
			lb::mmap_handle <char> ref_handle;
			ref_handle.open(reference_path);
			
			// Read the input FASTA.
			v2m::read_single_fasta_seq(ref_handle, reference, reference_seq_name);
			
			// Read the intermediate graph.
			cereal::PortableBinaryInputArchive iarchive(input_graph_stream);
			iarchive(graph);
		}
		
		// Output.
		auto const &aligned_ref_positions(graph.aligned_ref_positions());
		auto const &subgraph_start_positions(graph.subgraph_start_positions());
		auto const &alt_edge_count_csum(graph.alt_edge_count_csum());
		auto const &alt_edge_targets(graph.alt_edge_targets());
		auto const &alt_edge_labels(graph.alt_edge_labels());
		output_graph_stream << "digraph variants {\n";
		output_graph_stream << "\trankdir = LR;\n";
		output_graph_stream << "\trank = same;\n";
		
		// REF edges.
		for (auto const &pair : aligned_ref_positions | ranges::view::drop(1) | ranges::view::sliding(2))
		{
			auto const lhs(pair[0]);
			auto const rhs(pair[1]);
			std::string_view const ref_sub(reference.data() + lhs, rhs - lhs);
			output_graph_stream << "\t" << lhs << " -> " << rhs << " [label = \"" << ref_sub << "\", penwidth = 2.0];\n"; // FIXME: handle special characters?
		}
		output_graph_stream << '\n';
		
		// ALT edges.
		{
			std::size_t i(0);
			for (auto const &pair : alt_edge_count_csum | ranges::view::sliding(2))
			{
				auto const start(pair[0]);
				auto const end(pair[1]);
				for (auto const &[target, label] : ranges::view::slice(ranges::view::zip(alt_edge_targets, alt_edge_labels), start, end))
				{
					auto const src(aligned_ref_positions[1 + i]);
					auto const dst(aligned_ref_positions[1 + target]);
					output_graph_stream << "\t" << src << " -> " << dst << " [label = \"" << label << "\"];\n"; // FIXME: handle special characters?
				}
				
				++i;
			}
		}
		output_graph_stream << "}\n";
	}
}


int main(int argc, char **argv)
{
#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif

	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	try
	{
		output_for_gv(
			args_info.reference_arg,
			args_info.variants_arg,
			args_info.output_arg,
			args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
			args_info.overwrite_flag
		);
	}
	catch (lb::assertion_failure_exception const &exc)
	{
		std::cerr << "Assertion failure: " << exc.what() << '\n';
		boost::stacktrace::stacktrace const *st(boost::get_error_info <lb::traced>(exc));
		if (st)
			std::cerr << "Stack trace:\n" << *st << '\n';
		throw exc;
	}
	
	return EXIT_SUCCESS;
}
