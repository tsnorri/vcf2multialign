/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <list>
#include <range/v3/view/sliding.hpp>
#include <unistd.h>
#include <vcf2multialign/preprocess/variant_graph.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef std::vector <lb::file_ostream> output_stream_vector;
	
	
	struct stream_position
	{
		std::size_t node{};
		std::size_t stream_number{};
		
		bool operator<(stream_position const &other) const { return node < other.node; }
	};
	
	
	void output_chunk(std::string_view const &reference_sv, v2m::variant_graph const &graph, output_stream_vector &output_files);
	
	
	void output_founders(
		char const *reference_path,
		char const *input_graph_path,
		char const *reference_seq_name,
		std::size_t const chunk_size,
		bool const may_overwrite
	)
	{
		v2m::vector_type reference;
		v2m::variant_graph graph;
		
		// Open the input files.
		{
			lb::file_istream input_graph_stream;
			lb::open_file_for_reading(input_graph_path, input_graph_stream);
			
			lb::mmap_handle <char> ref_handle;
			ref_handle.open(reference_path);
			
			// Read the input FASTA.
			v2m::read_single_fasta_seq(ref_handle, reference, reference_seq_name);
			
			// Read the intermediate graph.
			cereal::PortableBinaryInputArchive iarchive(input_graph_stream);
			iarchive(graph);
		}
		
		// Create a string view from the reference.
		std::string_view const reference_sv(reference.data(), reference.size());
		
		// Output in chunks.
		{
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(may_overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			
			auto const max_paths(graph.max_paths_in_subgraph());
			auto const stream_count(1 + max_paths);
			std::size_t const chunk_count(std::ceil(1.0 * stream_count / chunk_size));
			for (auto const &pair : ranges::view::closed_iota(std::size_t(0), chunk_count) | ranges::view::sliding(2))
			{
				auto const lhs(pair[0]);
				auto const rhs(pair[1]);
				auto const rhs_(std::min(rhs, stream_count));
				std::cerr << "Processing chunk " << rhs << '/' << chunk_count << "…\n";
				output_stream_vector output_files(chunk_size * (rhs_ - lhs)); // Cannot reuse b.c. lb::file_ostream has a deleted copy constructor.
				
				// Open the output files.
				// Handle REF.
				if (0 == lhs)
				{
					lb::open_file_for_writing("REF", output_files.front(), mode);
					for (auto &&[i, of] : ranges::view::zip(ranges::view::iota(chunk_size * lhs, chunk_size * rhs_), output_files) | ranges::view::tail)
						lb::open_file_for_writing(std::to_string(i), of, mode);
				}
				else
				{
					for (auto &&[i, of] : ranges::view::zip(ranges::view::iota(chunk_size * lhs, chunk_size * rhs_), output_files))
						lb::open_file_for_writing(std::to_string(i), of, mode);
				}
				
				output_chunk(reference_sv, graph, output_files);
			}
		}
	}
	
	
	void output_chunk(std::string_view const &reference_sv, v2m::variant_graph const &graph, output_stream_vector &output_files)
	{
		// Output.
		{
			typedef v2m::variant_graph::sample_path_vector	sample_paths_type;
			typedef v2m::variant_graph::path_edge_matrix	path_edges_type;
			
			auto const &ref_positions(graph.ref_positions());
			auto const &aln_positions(graph.aligned_ref_positions());
			auto const &subgraph_start_positions(graph.subgraph_start_positions());
			auto const &alt_edge_count_csum(graph.alt_edge_count_csum());
			auto const &alt_edge_targets(graph.alt_edge_targets());
			auto const &alt_edge_labels(graph.alt_edge_labels());
			auto const &sample_paths(graph.sample_paths());
			auto const &path_edges(graph.path_edges());
			
			// Use a simple queue for keeping track of the aligned positions of the founders.
			// first: node index, second: stream number.
			std::list <stream_position> files_available(output_files.size()), files_waiting;
			for (auto &&[i, sp] : ranges::view::enumerate(files_available))
				sp.stream_number = i;
			
			// Handle the subgraphs, include the final subgraph.
			// Using ranges::view::single instead of repeat_n causes errors when compiling.
			// FIXME: normalize the file format s.t. the special cases in rsv are not needed.
			bool const first_subgraph_starts_from_zero(0 == subgraph_start_positions.front());
			auto const rsv(
				ranges::view::enumerate(
					ranges::view::zip(
						ranges::view::concat(
							ranges::view::repeat_n(sample_paths_type(), (first_subgraph_starts_from_zero ? 0 : 1)),
							sample_paths,
							ranges::view::repeat(sample_paths_type())
						),
						ranges::view::concat(
							ranges::view::repeat_n(path_edges_type(), (first_subgraph_starts_from_zero ? 0 : 1)),
							path_edges,
							ranges::view::repeat(path_edges_type())
						),
						ranges::view::concat(
							ranges::view::repeat_n(0, (first_subgraph_starts_from_zero ? 0 : 1)),	// Handle the possible padding. (The initial zero could be included when creating the graph.)
							subgraph_start_positions,
							ranges::view::single(ref_positions.size() - 2)
						) | ranges::view::sliding(2)
					)
				)
			);
			for (auto const &[subgraph_idx, tup] : rsv)
			{
				std::cerr << "Subgraph " << (1 + subgraph_idx) << '/' << (subgraph_start_positions.size() + (first_subgraph_starts_from_zero ? 0 : 1)) << "…\n";
				
				auto const &[sample_paths, edges_by_path_and_variant, ssp_pair] = tup;
				auto const subgraph_lhs(ssp_pair[0]);
				auto const subgraph_rhs(ssp_pair[1]);
				auto const subgraph_variants(edges_by_path_and_variant.number_of_rows());
				auto const subgraph_paths(edges_by_path_and_variant.number_of_columns());
				std::size_t subgraph_variant_idx(0); // within the subgraph.
				for (auto const &[node_idx, idx_pair] : ranges::view::enumerate(ranges::view::closed_iota(subgraph_lhs, subgraph_rhs) | ranges::view::sliding(2)))
				{
					std::cerr << "Node " << (1 + node_idx) << '/' << (subgraph_rhs - subgraph_lhs) << "…\n";
					
					auto const lhs(idx_pair[0]);
					auto const rhs(idx_pair[1]);
					auto const alt_lhs(alt_edge_count_csum[lhs]);
					auto const alt_rhs(alt_edge_count_csum[rhs]);
					auto const ref_lhs(ref_positions[1 + lhs]);
					auto const ref_rhs(ref_positions[1 + rhs]);
					auto const aln_lhs(aln_positions[1 + lhs]);
					auto const aln_rhs(aln_positions[1 + rhs]);
					auto const ref_len(ref_rhs - ref_lhs);
					auto const aln_len(aln_rhs - aln_lhs);
					
					// Find the range of output streams the current output position of which is before the current REF node index.
					// We iterate over the found files anyway, so no need to attempt sub-linear time.
					auto const it(std::find_if(files_available.begin(), files_available.end(), [lhs](auto const &sp){ return lhs < sp.node; }));
					
					// Move the items to the working list.
					libbio_assert(files_waiting.empty());
					files_waiting.splice(files_waiting.end(), files_available, files_available.begin(), it);
					
					// Write the sequences.
					if (0 == alt_rhs - alt_lhs)
					{
						// No ALT edges.
						// It may be that ref_len ≠ aln_len, though, if the rhs node has an in-ALT-edge.
						for (auto &sp : files_waiting)
						{
							auto const ref_begin(ref_positions[1 + sp.node]); // Start from the current node position.
							auto const ref_sub(reference_sv.substr(ref_begin, ref_rhs - ref_begin));
							auto &stream(output_files[sp.stream_number]);
							auto const t1(stream.tellp());
							stream << ref_sub;
							std::fill_n(std::ostream_iterator <char>(stream), aln_len - ref_len, '-');
							auto const t2(stream.tellp());
							libbio_assert_eq(aln_rhs, stream.tellp());
							sp.node = rhs;
						}
					}
					else
					{
						// Associate the path numbers with stream numbers directly (i.e. no weight-based matching).
						for (auto &sp : files_waiting)
						{
							auto &stream(output_files[sp.stream_number]);
							std::size_t const edge_number(sp.stream_number < subgraph_paths ? edges_by_path_and_variant(subgraph_variant_idx, sp.stream_number) : 0);
							if (0 == edge_number)
							{
								// REF edge.
								auto const ref_begin(ref_positions[1 + sp.node]);
								auto const ref_sub(reference_sv.substr(ref_lhs, ref_len));
								stream << ref_sub;
								std::fill_n(std::ostream_iterator <char>(stream), aln_len - ref_len, '-');
								libbio_assert_eq(aln_rhs, stream.tellp());
								sp.node = rhs;
							}
							else
							{
								// ALT edge.
								auto const alt_idx(alt_lhs + edge_number - 1);
								auto const target_node(alt_edge_targets[alt_idx]);
								auto const &alt(alt_edge_labels[alt_idx]);
								auto const aln_lhs(aln_positions[1 + lhs]);
								auto const aln_rhs(aln_positions[1 + target_node]);
								auto const aln_len(aln_rhs - aln_lhs);
								stream << alt;
								std::fill_n(std::ostream_iterator <char>(stream), aln_len - alt.size(), '-');
								libbio_assert_eq(aln_rhs, stream.tellp());
								sp.node = target_node;
							}
						}

						++subgraph_variant_idx;
					}
					
					// Sort the stream numbers by the writing position, i.e. node number.
					files_waiting.sort();
					files_available.merge(files_waiting);
				}
			}
			
			// Fill with the reference up to aligned reference length.
			{
				std::cerr << "Filling with the reference…\n";
				auto const ref_end(ref_positions.back());
				for (auto const &sp : files_available)
				{
					// Fill.
					auto const ref_begin(ref_positions[1 + sp.node]);
					auto const ref_sub(reference_sv.substr(ref_begin, ref_end - ref_begin));
					auto &stream(output_files[sp.stream_number]);
					stream << ref_sub << std::flush;
				}
			}
			
			std::cerr << "Finishing…\n";
		}
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
	
	if (args_info.chunk_size_arg <= 0)
	{
		std::cerr << "Chunk size must be positive." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	try
	{
		output_founders(
			args_info.reference_arg,
			args_info.variants_arg,
			args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
			args_info.chunk_size_arg,
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
