/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/subprocess.hh>
#include <map>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <string>
#include <string_view>
#include <vcf2multialign/variant_graph.hh>
#include <vector>
#include "cmdline.h"

namespace ios	= boost::iostreams;
namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef lb::subprocess <lb::subprocess_handle_spec::STDIN>	subprocess_type;
	
	
	void open_stream_with_file_handle(
		lb::file_ostream &stream,
		lb::file_handle &fh
	)
	{
		stream.open(fh.get(), ios::never_close_handle);
		stream.exceptions(std::ostream::badbit);
	}
	
	
	void output_graphviz(
		v2m::sequence_type const &ref_seq_,
		v2m::variant_graph const &graph,
		lb::file_handle &fh
	)
	{
		std::string_view const ref_seq(ref_seq_.data(), ref_seq_.size());
		
		lb::file_ostream stream;
		stream.open(fh.get(), ios::never_close_handle);
		stream.exceptions(std::ostream::badbit);
		
		stream << "digraph variants {\n";
		stream << "\trankdir = LR;\n";
		stream << "\trank = same;\n";
		
		typedef v2m::variant_graph				variant_graph;
		typedef variant_graph::position_type	position_type;
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::edge_type		edge_type;
		
		// Nodes.
		for (auto const &[node, ref_pos, aln_pos] : rsv::zip(rsv::iota(0), graph.reference_positions, graph.aligned_positions))
			stream << '\t' << node << " [shape = Mrecord, label = \"" << node << " | " << ref_pos << " | " << aln_pos << "\"];\n";
		stream << '\n';
		
		// REF edges.
		for (auto const &[node, range] : rsv::enumerate(graph.reference_positions | rsv::sliding(2)))
		{
			auto const lb(range[0]);
			auto const rb(range[1]);
			auto const label(ref_seq.substr(lb, rb - lb));

			// FIXME: Handle special characters in the label?
			stream << '\t' << node << " -> " << (node + 1) << " [label = \"";
			if (label.size() <= label.size())
				stream << label;
			else
				stream << label.substr(0, 10) << "…" << label.substr(label.size() - 10, 10) << " (" << label.size() << ')';
			stream << "\", penwidth = 2.0];\n";
		}
		stream << '\n';
		
		// ALT edges.
		for (auto const &[src_node, edge_range] : rsv::enumerate(graph.alt_edge_count_csum | rsv::sliding(2)))
		{
			auto const edge_lb(edge_range[0]);
			auto const edge_rb(edge_range[1]);
			
			for (edge_type edge_idx(edge_lb); edge_idx < edge_rb; ++edge_idx)
				stream << '\t' << src_node << " -> " << graph.alt_edge_targets[edge_idx] << " [label = \"" << graph.alt_edge_labels[edge_idx] << "\"];\n";
		}
		stream << "}\n";
	}
	
	
	void output_sequence(
		v2m::sequence_type const &ref_seq,
		v2m::variant_graph const &graph,
		v2m::variant_graph::sample_type const sample_idx,
		v2m::variant_graph::ploidy_type const chr_copy_idx,
		lb::file_ostream &stream
	)
	{
		typedef v2m::variant_graph				variant_graph;
		typedef variant_graph::position_type	position_type;
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::edge_type		edge_type;
		
		position_type ref_pos{};
		position_type aln_pos{};
		position_type next_ref_pos{};
		position_type next_aln_pos{};
		node_type current_node{};
		auto const limit(graph.node_count() - 1);
		auto const chr_copy_idx_(variant_graph::SAMPLE_MAX == sample_idx ? 0 : graph.ploidy_csum[sample_idx] + chr_copy_idx);
		while (current_node < limit)
		{
			std::size_t label_size{};
			if (variant_graph::SAMPLE_MAX != sample_idx) // Always follow REF edges if outputting the aligned reference.
			{
				auto const &[edge_lb, edge_rb] = graph.edge_range_for_node(current_node);
				for (edge_type edge_idx(edge_lb); edge_idx < edge_rb; ++edge_idx)
				{
					if (graph.paths_by_chrom_copy_and_edge(edge_idx, chr_copy_idx_))
					{
						// Found an ALT edge to follow.
						auto const target_node(graph.alt_edge_targets[edge_idx]);
						auto const &label(graph.alt_edge_labels[edge_idx]);
						next_ref_pos = graph.reference_positions[target_node];
						next_aln_pos = graph.aligned_positions[target_node];
						libbio_assert_lte(label.size(), next_aln_pos - aln_pos);
						stream << label;
						current_node = target_node;
						label_size = label.size();
						goto continue_loop;
					}
				}
			}
			
			{
				next_ref_pos = graph.reference_positions[current_node + 1];
				next_aln_pos = graph.aligned_positions[current_node + 1];
				std::string_view const ref_part(ref_seq.data() + ref_pos, next_ref_pos - ref_pos);
				stream << ref_part;
				label_size = ref_part.size();
				++current_node;
			}
			
		continue_loop:
			std::fill_n(std::ostreambuf_iterator <char>(stream), next_aln_pos - aln_pos - label_size, '-');
			ref_pos = next_ref_pos;
			aln_pos = next_aln_pos;
		}
	}
	
	
	void output_sequence(
		v2m::sequence_type const &ref_seq,
		v2m::variant_graph const &graph,
		v2m::variant_graph::sample_type const sample_idx,
		v2m::variant_graph::ploidy_type const chr_copy_idx,
		lb::file_handle &fh
	)
	{
		lb::file_ostream stream;
		open_stream_with_file_handle(stream, fh);
		output_sequence(ref_seq, graph, sample_idx, chr_copy_idx, stream);
	}
	
	
	void exit_subprocess(subprocess_type &proc)
	{
		auto const res(proc.close());
		auto const &[close_status, exit_status, pid] = res;
		if (! (lb::process_handle::close_status::exit_called == close_status && 0 == exit_status))
		{
			// Try to determine the reason for the exit status.
			auto const &status(proc.status());
			switch (status.execution_status)
			{
				case lb::execution_status_type::no_error:
				{
					std::cerr << "ERROR: Subprocess with PID " << pid << " exited with status " << exit_status;
					switch (close_status)
					{
						case lb::process_handle::close_status::unknown:
							std::cerr << " (exiting reason not known)";
							break;
						case lb::process_handle::close_status::terminated_by_signal:
							std::cerr << " (terminated by signal)";
							break;
						case lb::process_handle::close_status::stopped_by_signal:
							std::cerr << " (stopped by signal)";
							break;
						default:
							break;
					}
					break;
				}
				
				case lb::execution_status_type::file_descriptor_handling_failed:
				case lb::execution_status_type::exec_failed:
				{
					std::cerr << "ERROR: Unable to start subprocess: " << strerror(status.error);
					break;
				}
			}
			
			std::cerr << '\n';
			std::exit(EXIT_FAILURE);
		}
	}
	
	
	void output_sequence_file(
		v2m::sequence_type const &ref_seq,
		v2m::variant_graph const &graph,
		v2m::variant_graph::sample_type const sample_idx,
		v2m::variant_graph::ploidy_type const chr_copy_idx,
		char const * const pipe_cmd,
		char const * const dst_name
	)
	{
		if (pipe_cmd)
		{
			auto proc(subprocess_type::subprocess_with_arguments({pipe_cmd, dst_name}));
			auto &fh(proc.stdin_handle());
			output_sequence(ref_seq, graph, sample_idx, chr_copy_idx, fh);
			exit_subprocess(proc);
		}
		else
		{
			lb::file_handle fh(lb::open_file_for_writing(dst_name, lb::writing_open_mode::CREATE));
			output_sequence(ref_seq, graph, sample_idx, chr_copy_idx, fh);
		}
	}
	
	
	void output_sequence_files(v2m::sequence_type const &ref_seq, v2m::variant_graph const &graph, char const * const pipe_cmd)
	{
		typedef v2m::variant_graph			variant_graph;
		typedef variant_graph::ploidy_type	ploidy_type;
		
		output_sequence_file(ref_seq, graph, variant_graph::SAMPLE_MAX, 0, pipe_cmd, "REF");
		for (auto const &[sample_idx, sample] : rsv::enumerate(graph.sample_names))
		{
			auto const ploidy(graph.sample_ploidy(sample_idx));
			for (auto const chr_copy_idx : rsv::iota(ploidy_type(0), ploidy))
			{
				// FIXME: Use std::format.
				std::stringstream dst_name;
				dst_name << sample;
				dst_name << '-';
				dst_name << chr_copy_idx;
				output_sequence_file(ref_seq, graph, sample_idx, chr_copy_idx, pipe_cmd, dst_name.str().data());
			}
		}
	}
	
	
	void output_sequences_a2m(v2m::sequence_type const &ref_seq, v2m::variant_graph const &graph, lb::file_ostream &stream)
	{
		typedef v2m::variant_graph			variant_graph;
		typedef variant_graph::ploidy_type	ploidy_type;
		
		stream << ">REF\n";
		output_sequence(ref_seq, graph, variant_graph::SAMPLE_MAX, 0, stream);
		stream << '\n';
		
		std::uint32_t seq_count{1};
		for (auto const &[sample_idx, sample] : rsv::enumerate(graph.sample_names))
		{
			auto const ploidy(graph.sample_ploidy(sample_idx));
			for (auto const chr_copy_idx : rsv::iota(ploidy_type(0), ploidy))
			{
				stream << '>' << sample << '-' << chr_copy_idx << '\n';
				output_sequence(ref_seq, graph, sample_idx, chr_copy_idx, stream);
				stream << '\n';

				++seq_count;
				if (0 == seq_count % 10)
					lb::log_time(std::cerr) << "Handled " << seq_count << " sequences…\n";
			}
		}
	}
	
	
	void output_sequences_a2m(v2m::sequence_type const &ref_seq, v2m::variant_graph const &graph, char const * const dst_name, char const * const pipe_cmd)
	{
		if (pipe_cmd)
		{
			auto proc(subprocess_type::subprocess_with_arguments({pipe_cmd, dst_name}));
			auto &fh(proc.stdin_handle());
			
			{
				lb::file_ostream stream;
				open_stream_with_file_handle(stream, fh);
				output_sequences_a2m(ref_seq, graph, stream);
			}
			
			exit_subprocess(proc);
		}
		else
		{
			lb::file_handle fh(lb::open_file_for_writing(dst_name, lb::writing_open_mode::CREATE));
			lb::file_ostream stream;
			open_stream_with_file_handle(stream, fh);
			output_sequences_a2m(ref_seq, graph, stream);
		}
	}
	
	
	void run(
		char const *reference_path,
		char const *variants_path,
		char const *ref_seq_id,
		char const *chr_id,
		char const *sequence_a2m_output_path,
		bool const should_output_sequences_separate,
		char const *pipe_cmd,
		char const *graphviz_output_path
	)
	{
		// Read the reference sequence.
		v2m::sequence_type ref_seq;
		{
			if (ref_seq_id)
				lb::log_time(std::cerr) << "Reading reference sequence with identifier “" << ref_seq_id << "”…" << std::flush;
			else
				lb::log_time(std::cerr) << "Reading the first reference sequence from the input FASTA…" << std::flush;
			auto const res(lb::read_single_fasta_sequence(reference_path, ref_seq, ref_seq_id));
			
			if (!res)
			{
				std::cerr << " ERROR: Unable to read the reference sequence.\n";
				std::exit(EXIT_FAILURE);
			}
			
			std::cerr << " Done. Reference length is " << ref_seq.size() << ".\n";
		}
		
		lb::log_time(std::cerr) << "Building the variant graph…\n";
		v2m::variant_graph graph;
		v2m::build_graph_statistics stats;
		v2m::build_variant_graph(ref_seq, variants_path, chr_id, graph, stats, &std::cout);
		lb::log_time(std::cerr) << "Done. Handled variants: " << stats.handled_variants << " chromosome ID mismatches: " << stats.chr_id_mismatches << "\n";
		
		if (graphviz_output_path)
		{
			lb::log_time(std::cerr) << "Outputting the variant graph in Grapnviz format…" << std::flush;
			lb::file_handle fh(lb::open_file_for_writing(graphviz_output_path, lb::writing_open_mode::CREATE));
			output_graphviz(ref_seq, graph, fh);
			std::cerr << " Done.\n";
		}
		
		if (sequence_a2m_output_path)
		{
			lb::log_time(std::cerr) << "Outputting sequences as A2M…\n";
			output_sequences_a2m(ref_seq, graph, sequence_a2m_output_path, pipe_cmd);
			lb::log_time(std::cerr) << "Done.\n";
		}
		
		if (should_output_sequences_separate)
		{
			lb::log_time(std::cerr) << "Outputting sequences one by one…" << std::flush;
			output_sequence_files(ref_seq, graph, pipe_cmd);
			std::cerr << " Done.\n";
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
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	if (args_info.show_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}
	
	try
	{
		run(
			args_info.reference_arg,
			args_info.variants_arg,
			args_info.reference_sequence_arg,
			args_info.chromosome_arg,
			args_info.output_sequences_a2m_arg,
			args_info.output_sequences_separate_flag,
			args_info.pipe_arg,
			args_info.output_graphviz_arg
		);
	}
	catch (std::exception const &exc)
	{
		std::cerr << "ERROR: Caught an exception: " << exc.what() << '\n';
		auto const trace(boost::get_error_info <lb::traced>(exc));
		if (trace)
			std::cerr << "Stack trace:\n" << (*trace) << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
