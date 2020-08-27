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
#include <range/v3/view/zip.hpp>
#include <unistd.h>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	void output_for_gv(
		char const *reference_path,
		char const *input_graph_path,
		char const *output_graph_path,
		char const *reference_seq_name,
		std::size_t const subgraph_idx,
		std::size_t max_aligned_pos,
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
		auto const &ref_positions(graph.ref_positions());
		auto const &aligned_ref_positions(graph.aligned_ref_positions());
		auto const &subgraph_start_positions(graph.subgraph_start_positions());
		auto const &alt_edge_count_csum(graph.alt_edge_count_csum());
		auto const &alt_edge_targets(graph.alt_edge_targets());
		auto const &alt_edge_labels(graph.alt_edge_labels());

		if (SIZE_MAX != subgraph_idx && subgraph_start_positions.size() <= subgraph_idx)
		{
			std::cerr << "Got subgraph index " << subgraph_idx << " while there are " << subgraph_start_positions.size() << " subgraphs.\n";
			std::exit(EXIT_FAILURE);
		}

		// Check whether iteration should be started from a subgraph.
		// Convert the indices to 1-based.
		auto const pos_slice_start(1 + (SIZE_MAX == subgraph_idx ? 0 : subgraph_start_positions[subgraph_idx]));
		auto const pos_slice_end(
			SIZE_MAX == subgraph_idx || subgraph_start_positions.size() <= (1 + subgraph_idx)
			? ref_positions.size()
			: 2 + subgraph_start_positions[1 + subgraph_idx]
		);

		// Update max_aligned_pos to point to the subgraph end.
		max_aligned_pos = std::min(max_aligned_pos, aligned_ref_positions[pos_slice_end - 1]);

		output_graph_stream << "digraph variants {\n";
		output_graph_stream << "\trankdir = LR;\n";
		output_graph_stream << "\trank = same;\n";
		
		// Labels if needed, otherwise use node ID as the label.
		{
			// Since background colours are needed at subgraph start positions, handle them, too.
			auto sg_it(subgraph_start_positions.begin() + (SIZE_MAX == subgraph_idx ? 0 : subgraph_idx));
			auto const sg_end(subgraph_start_positions.end());
			std::size_t next_subgraph_start_idx(sg_it == sg_end ? SIZE_MAX : *sg_it++);

			auto const range(
				rsv::zip(
					rsv::ints(0),
					ref_positions | rsv::slice(pos_slice_start, pos_slice_end),
					aligned_ref_positions | rsv::slice(pos_slice_start, pos_slice_end)
				)
			);
			for (auto const &[i, ref_pos, aligned_ref_pos] : range)
			{
				if (max_aligned_pos < aligned_ref_pos)
					break;

				if (i == next_subgraph_start_idx)
				{
					next_subgraph_start_idx = (sg_it == sg_end ? SIZE_MAX : *sg_it++);
					output_graph_stream << '\t' << aligned_ref_pos << " [shape = Mrecord, label = \"" << ref_pos << " | " << aligned_ref_pos << "\", style = filled, fillcolor = grey95];\n";
				}
				else
				{
					output_graph_stream << '\t' << aligned_ref_pos << " [shape = Mrecord, label = \"" << ref_pos << " | " << aligned_ref_pos << "\"];\n";
				}
			}
		}
		output_graph_stream << '\n';

		// REF edges.
		{
			auto const range(rsv::zip(
				ref_positions | rsv::slice(pos_slice_start, pos_slice_end) | rsv::sliding(2),
				aligned_ref_positions | rsv::slice(pos_slice_start, pos_slice_end) | rsv::sliding(2)
			));
			for (auto const &[ref_pair, aligned_ref_pair] : range)
			{
				auto const aligned_lhs(aligned_ref_pair[0]);
				auto const aligned_rhs(aligned_ref_pair[1]);

				if (max_aligned_pos < aligned_lhs)
					break;

				if (max_aligned_pos < aligned_rhs)
					continue;

				auto const ref_lhs(ref_pair[0]);
				auto const ref_rhs(ref_pair[1]);
				std::string_view const ref_sub(reference.data() + ref_lhs, ref_rhs - ref_lhs);
				
				output_graph_stream << '\t' << aligned_lhs << " -> " << aligned_rhs << " [label = \"";
				// If the text is longer than 20 characters, truncate and add an ellipsis.
				// FIXME: handle special characters?
				if (20 < ref_sub.size())
				{
					output_graph_stream << ref_sub.substr(0, 10);
					output_graph_stream << "…";
					output_graph_stream << ref_sub.substr(ref_sub.size() - 10, 10);
					output_graph_stream << " (" << ref_sub.size() << ')';
				}
				else
				{
					output_graph_stream << ref_sub;
				}
				output_graph_stream << "\", penwidth = 2.0];\n";
			}
			output_graph_stream << '\n';
		}
		
		// ALT edges.
		{
			std::size_t i(pos_slice_start);
			for (auto const &pair : alt_edge_count_csum | rsv::slice(pos_slice_start, pos_slice_end) | rsv::sliding(2))
			{
				auto const start(pair[0]);
				auto const end(pair[1]);
				for (auto const &[target, label] : rsv::slice(rsv::zip(alt_edge_targets, alt_edge_labels), start, end))
				{
					auto const src(aligned_ref_positions[1 + i]);
					auto const dst(aligned_ref_positions[1 + target]);

					if (max_aligned_pos < src)
						break;

					if (max_aligned_pos < dst)
						continue;

					output_graph_stream << '\t' << src << " -> " << dst << " [label = \"" << label << "\"];\n"; // FIXME: handle special characters?
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
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

	std::size_t subgraph_idx(SIZE_MAX);
	if (args_info.subgraph_given)
	{
		if (args_info.subgraph_arg < 0)
		{
			std::cerr << "Subgraph index needs to be non-negative.\n";
			std::exit(EXIT_FAILURE);
		}
		subgraph_idx = args_info.subgraph_arg;
	}

	std::size_t max_aligned_pos(SIZE_MAX);
	if (args_info.max_aligned_pos_given)
	{
		if (args_info.max_aligned_pos_arg < 0)
		{
			std::cerr << "Maximum aligned position needs to be non-negative.\n";
			std::exit(EXIT_FAILURE);
		}
		max_aligned_pos = args_info.max_aligned_pos_arg;
	}
	
	try
	{
		output_for_gv(
			args_info.reference_arg,
			args_info.variants_arg,
			args_info.output_arg,
			args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
			subgraph_idx,
			max_aligned_pos,
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
