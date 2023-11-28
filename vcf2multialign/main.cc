/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/generic_parser.hh>
#include <libbio/subprocess.hh>
#include <map>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <string>
#include <string_view>
#include <vcf2multialign/output.hh>
#include <vcf2multialign/variant_graph.hh>
#include <vector>
#include "cmdline.h"

namespace ios	= boost::iostreams;
namespace lb	= libbio;
namespace lbp	= libbio::parsing;
namespace rsv	= ranges::views;
namespace v2m	= vcf2multialign;


namespace {
	
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
	
	
	template <typename t_string>
	struct sample_identifier_tpl
	{
		t_string		sample;
		std::uint32_t	chromosome_copy_index{};

		template <typename t_other_string>
		sample_identifier_tpl(t_other_string &&sample_, std::uint32_t const chromosome_copy_index_):
			sample(std::forward <t_other_string>(sample_)),
			chromosome_copy_index(chromosome_copy_index_)
		{
		}

		decltype(auto) to_tuple() const { return std::tie(sample, chromosome_copy_index); }
		bool operator<(sample_identifier_tpl const &other) const { return to_tuple() < other.to_tuple(); }
		
		template <typename t_other_string>
		bool operator<(sample_identifier_tpl <t_other_string> const &other) const { return to_tuple() < other.to_tuple(); }
	};

	typedef sample_identifier_tpl <std::string>			sample_identifier;
	typedef sample_identifier_tpl <std::string_view>	sample_identifier_sv;


	struct build_variant_graph_delegate final : public v2m::build_graph_delegate
	{
		std::vector <sample_identifier>	excluded_samples;

		void report_overlapping_alternative(
			std::string_view const sample_name,
			v2m::variant_graph::ploidy_type const chrom_copy_idx,
			v2m::variant_graph::position_type const ref_pos,
			std::vector <std::string_view> const &var_id,
			std::uint32_t const gt
		) override
		{
			std::cout << "Overlapping alternative alleles. Sample: " << sample_name << " chromosome copy: " << chrom_copy_idx << " current variant position: " << ref_pos << " genotype: " << gt << '\n';
		}

		bool should_include(std::string_view const sample_name, v2m::variant_graph::ploidy_type const chrom_copy_idx) const override
		{
			return !std::binary_search(excluded_samples.begin(), excluded_samples.end(), sample_identifier_sv{sample_name, chrom_copy_idx});
		}
		
		bool ref_column_mismatch(std::uint64_t const var_idx, position_type const pos, std::string_view const expected, std::string_view const actual) override
		{
			std::cerr << "WARNING: REF column contents do not match the reference sequence in variant " << var_idx << ". Expected: “" << expected << "” Actual: “" << actual << "”\n";
			return true;
		}
	};


	void read_excluded_samples(char const *exclude_samples_tsv_path, char const *chr_id, std::vector <sample_identifier> &excluded_samples)
	{
		lb::file_istream stream;
		lb::open_file_for_reading(exclude_samples_tsv_path, stream);

		typedef lbp::parser <
			lbp::traits::delimited <lbp::delimiter <'\t'>, lbp::delimiter <'\n'>>,
			lbp::fields::text <>,
			lbp::fields::text <>,
			lbp::fields::integer <std::uint32_t>
		> parser_type;
		typedef parser_type::record_type record_type;

		{
			std::istreambuf_iterator it(stream.rdbuf());
			std::istreambuf_iterator <decltype(it)::value_type> const sentinel;
			parser_type parser;
			record_type rec;
			try
			{
				while (true)
				{
					if (!parser.parse(it, sentinel, rec))
						break;

					if (chr_id == std::get <0>(rec))
						excluded_samples.emplace_back(std::get <1>(rec), std::get <2>(rec));
				}
			}
			catch (lbp::parse_error const &err)
			{
				std::cerr << "ERROR: Parse error at position " << err.position() << ": ";
				err.output_error(std::cerr);
				std::cerr << '\n';
				std::exit(EXIT_FAILURE);
			}
		}

		std::sort(excluded_samples.begin(), excluded_samples.end());
	}
	
	
	void build_variant_graph(
		char const *variants_path,
		char const *chr_id,
		char const *exclude_samples_tsv_path,
		v2m::sequence_type const &ref_seq,
		v2m::variant_graph &graph,
		bool const be_verbose
	)
	{
		build_variant_graph_delegate delegate;
		if (exclude_samples_tsv_path)
		{
			lb::log_time(std::cerr) << "Reading the excluded sample list…" << std::flush;
			read_excluded_samples(exclude_samples_tsv_path, chr_id, delegate.excluded_samples);
			std::cerr << " Done.\n";
			
			if (be_verbose)
			{
				std::cerr << "Excluded the following samples:\n";
				for (auto const &sample_id : delegate.excluded_samples)
					std::cerr << sample_id.sample << " (" << sample_id.chromosome_copy_index << ")\n";
			}
		}
		
		lb::log_time(std::cerr) << "Building the variant graph…\n";
		v2m::build_graph_statistics stats;
		v2m::build_variant_graph(ref_seq, variants_path, chr_id, graph, stats, delegate);
		lb::log_time(std::cerr) << "Done. Handled variants: " << stats.handled_variants << " chromosome ID mismatches: " << stats.chr_id_mismatches << "\n";
	}
	
	
	class output_delegate final : public v2m::output_delegate
	{
	private:
		v2m::variant_graph const	*m_graph{};
		bool						m_is_verbose{};
		
	public:
		output_delegate(v2m::variant_graph const &graph, bool const is_verbose):
			m_graph(&graph),
			m_is_verbose(is_verbose)
		{
		}
		
		
		void will_handle_sample(std::string const &sample, sample_type const sample_idx, ploidy_type const chr_copy_idx) override
		{
			if (m_is_verbose)
				lb::log_time(std::cerr) << "Sample: " << sample << " (" << (1 + sample_idx) << "/" << m_graph->sample_names.size() << ") copy index: " << chr_copy_idx << '\n';
		}
		
		
		void will_handle_founder_sequence(sample_type const sample_idx) override
		{
			if (m_is_verbose)
				lb::log_time(std::cerr) << "Founder sequence " << sample_idx << '\n';
		}
		
		
		void handled_sequences(sequence_count_type const seq_count) override
		{
			if (0 == seq_count % 10)
			{
				auto const total_seq_count(m_graph->total_chromosome_copies());
				lb::log_time(std::cerr) << "Handled " << seq_count << '/' << total_seq_count << " sequences…\n";
			}
		}
		
		
		void exit_subprocess(v2m::subprocess_type &proc) override
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
	};
	
	
	void run(gengetopt_args_info const &args_info)
	{
		// Read the reference sequence.
		v2m::sequence_type ref_seq;
		{
			if (args_info.reference_sequence_arg)
				lb::log_time(std::cerr) << "Reading reference sequence with identifier “" << args_info.reference_sequence_arg << "”…" << std::flush;
			else
				lb::log_time(std::cerr) << "Reading the first reference sequence from the input FASTA…" << std::flush;
			auto const res(lb::read_single_fasta_sequence(args_info.input_reference_arg, ref_seq, args_info.reference_sequence_arg));
			
			if (!res)
			{
				std::cerr << " ERROR: Unable to read the reference sequence.\n";
				std::exit(EXIT_FAILURE);
			}
			
			std::cerr << " Done. Reference length is " << ref_seq.size() << ".\n";
		}
		
		v2m::variant_graph graph;
		if (args_info.input_graph_given)
		{
			lb::log_time(std::cerr) << "Loading the variant graph from " << args_info.input_graph_arg << "…" << std::flush;
			lb::file_istream is;
			lb::open_file_for_reading(args_info.input_graph_arg, is);
			cereal::PortableBinaryInputArchive archive(is);
			archive(graph);
			std::cerr << " Done.\n";
		}
		else
		{
			build_variant_graph(args_info.input_variants_arg, args_info.chromosome_arg, args_info.exclude_samples_arg, ref_seq, graph, args_info.verbose_given);
		}
		
		if (args_info.output_graph_given)
		{
			lb::log_time(std::cerr) << "Outputting the variant graph…" << std::flush;
			lb::file_ostream os;
			lb::open_file_for_writing(args_info.output_graph_arg, os, lb::writing_open_mode::CREATE);
			cereal::PortableBinaryOutputArchive archive(os);
			archive(graph);
			std::cerr << " Done.\n";
		}
		
		if (args_info.output_graphviz_given)
		{
			lb::log_time(std::cerr) << "Outputting the variant graph in Graphviz format…" << std::flush;
			lb::file_handle fh(lb::open_file_for_writing(args_info.output_graphviz_arg, lb::writing_open_mode::CREATE));
			output_graphviz(ref_seq, graph, fh);
			std::cerr << " Done.\n";
		}
		
		{
			output_delegate delegate(graph, args_info.verbose_given);
			auto do_output([&args_info, &ref_seq, &graph](v2m::output &output){
				if (args_info.output_sequences_a2m_given)
				{
					lb::log_time(std::cerr) << "Outputting sequences as A2M…\n";
					output.output_a2m(ref_seq, graph, args_info.output_sequences_a2m_arg);
					lb::log_time(std::cerr) << "Done.\n";
				}
				
				if (args_info.output_sequences_separate_given)
				{
					lb::log_time(std::cerr) << "Outputting sequences one by one…" << std::flush;
					output.output_separate(ref_seq, graph, separate_output_format_arg_A2M == args_info.separate_output_format_arg);
					std::cerr << " Done.\n";
				}
			});
			
			if (args_info.haplotypes_given)
			{
				v2m::haplotype_output output(args_info.pipe_arg, args_info.unaligned_output_given, delegate);
				do_output(output);
			}
			else if (args_info.founder_sequences_given)
			{
				v2m::founder_sequence_greedy_output output(args_info.pipe_arg, args_info.keep_ref_edges_given, args_info.unaligned_output_given, delegate);
				
				if (args_info.input_cut_positions_given)
					output.load_cut_positions(args_info.input_cut_positions_arg);
				else
				{
					lb::log_time(std::cerr) << "Optimising cut positions…\n";
					if (!output.find_cut_positions(graph, args_info.minimum_distance_arg))
					{
						std::cerr << "ERROR: Unable to optimise cut positions.\n";
						std::exit(EXIT_FAILURE);
					}
					
					if (args_info.verbose_flag)
					{
						std::cout << "Cut positions:";
						for (auto const cp : output.cut_positions())
							std::cout << ' ' << cp;
						std::cout << '\n';
					}
				}

				std::cout << "Maximum segmentation height: " << (1 + output.max_segmentation_height()) << '\n';
				
				if (args_info.output_cut_positions_given)
					output.output_cut_positions(args_info.output_cut_positions_arg);
				
				lb::log_time(std::cerr) << "Finding matchings in the variant graph…\n";
				if (!output.find_matchings(graph, args_info.founder_sequences_arg))
				{
					std::cerr << "ERROR: Unable to find matchings.\n";
					std::exit(EXIT_FAILURE);
				}
				
				if (args_info.verbose_flag)
				{
					std::cout << "Matchings:\n";
					auto const &assigned_samples(output.assigned_samples());
					for (auto const col_idx : rsv::iota(std::size_t(0), assigned_samples.number_of_columns()))
					{
						auto const col(assigned_samples.column(col_idx));
						std::cout << col_idx << ':';
						for (auto const val : col)
							std::cout << '\t' << val;
						std::cout << '\n';
					}
				}
				
				do_output(output);
			}
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
	
	if (args_info.show_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}
	
	if (args_info.input_variants_given && args_info.input_graph_given)
	{
		std::cerr << "ERROR: Only one of --input-variants and --input-graph can be specified.\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (! (args_info.input_variants_given || args_info.input_graph_given))
	{
		std::cerr << "ERROR: One of --input-variants and --input-graph must be specified.\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (args_info.input_variants_given && !args_info.chromosome_given)
	{
		std::cerr << "ERROR: --chromosome must be specified with --input-variants.\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (args_info.founder_sequences_given && args_info.founder_sequences_arg <= 0)
	{
		std::cerr << "ERROR: --founder-sequences must be positive.\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (args_info.minimum_distance_given && args_info.input_cut_positions_given)
	{
		std::cerr << "ERROR: --input-cut-positions and --minimum-distance are mutually exclusive.\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (args_info.minimum_distance_arg < 0)
	{
		std::cerr << "ERROR: --minimum-distance must be non-negative.\n";
		std::exit(EXIT_FAILURE);
	}
	
	try
	{
		run(args_info);
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
