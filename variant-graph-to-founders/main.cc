/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <unistd.h>
#include <vcf2multialign/preprocess/variant_graph.hh>
#include "cmdline.h"
#include "sequence_generator_base.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	class founder_alt_edge_handler final : public v2m::alt_edge_handler_base
	{
	public:
		typedef v2m::alt_edge_handler_base::stream_position_list	stream_position_list;
		
		using v2m::alt_edge_handler_base::alt_edge_handler_base;
		
		void handle_node_with_alt_edges(stream_position_list &files_waiting, std::size_t const subgraph_variant_idx) const override
		{
			// Associate the path numbers with stream numbers directly (i.e. no weight-based matching).
			auto const subgraph_paths(m_edges_by_path_and_variant->number_of_columns());
			for (auto &sp : files_waiting)
			{
				auto &stream((*m_output_files)[sp.stream_number]);
				std::size_t const edge_number(sp.stream_number < subgraph_paths ? (*m_edges_by_path_and_variant)(subgraph_variant_idx, sp.stream_number) : 0);
				this->handle_edge(sp, edge_number);
			}
		}
	};
	
	
	class haplotype_alt_edge_handler final : public v2m::alt_edge_handler_base
	{
	public:
		typedef v2m::alt_edge_handler_base::stream_position_list	stream_position_list;
		
		using v2m::alt_edge_handler_base::alt_edge_handler_base;
		
		void handle_node_with_alt_edges(stream_position_list &files_waiting, std::size_t const subgraph_variant_idx) const override
		{
			this->handle_edge(files_waiting.back(), 0); // REF
			for (auto &sp : files_waiting | ranges::view::drop_last(1))
			{
				auto &stream((*m_output_files)[sp.stream_number]);
				auto const path_number((*m_sample_paths)[sp.stream_number]);
				auto const edge_number((*m_edges_by_path_and_variant)(subgraph_variant_idx, path_number));
				this->handle_edge(sp, edge_number);
			}
		}
	};
	
	
	class founder_sequence_generator final : public v2m::sequence_generator_base
	{
	public:
		typedef v2m::sequence_generator_base::output_stream_vector	output_stream_vector;
		
		std::unique_ptr <v2m::alt_edge_handler_base> make_alt_edge_handler(
			std::string_view const &reference_sv,
			v2m::variant_graph const &graph,
			output_stream_vector &output_files
		) const override
		{
			return std::unique_ptr <v2m::alt_edge_handler_base>(new founder_alt_edge_handler(reference_sv, graph, output_files));
		}
		
		std::size_t const get_stream_count(v2m::variant_graph const &graph) const override
		{
			auto const max_paths(graph.max_paths_in_subgraph());
			return 1 + max_paths;
		}
		
	protected:
		void open_output_file(std::size_t const idx, lb::file_ostream &of, lb::writing_open_mode const mode, v2m::variant_graph const &graph) const override
		{
			lb::open_file_for_writing(std::to_string(1 + idx), of, mode);
		}
	};
	
	
	class haplotype_sequence_generator final : public v2m::sequence_generator_base
	{
	public:
		typedef v2m::sequence_generator_base::output_stream_vector	output_stream_vector;
		
		std::unique_ptr <v2m::alt_edge_handler_base> make_alt_edge_handler(
			std::string_view const &reference_sv,
			v2m::variant_graph const &graph,
			output_stream_vector &output_files
		) const override
		{
			return std::unique_ptr <v2m::alt_edge_handler_base>(new haplotype_alt_edge_handler(reference_sv, graph, output_files));
		}
		
		std::size_t const get_stream_count(v2m::variant_graph const &graph) const override
		{
			return 1 + graph.sample_names().size();
		}
		
	protected:
		void open_output_file(std::size_t const idx, lb::file_ostream &of, lb::writing_open_mode const mode, v2m::variant_graph const &graph) const override
		{
			auto const &name(graph.sample_names()[idx]);
			lb::open_file_for_writing(name, of, mode);
		}
	};
	
	
	void output_sequences(v2m::sequence_generator_base &gen, gengetopt_args_info const &args_info)
	{
		gen.output_sequences(
			args_info.reference_arg,
			args_info.variants_arg,
			args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
			args_info.chunk_size_arg,
			args_info.overwrite_flag
		);
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
	
	{
		auto const mode_count(args_info.output_founders_given + args_info.output_haplotypes_given);
		if (0 == mode_count)
		{
			std::cerr << "No mode given." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	try
	{
		if (args_info.output_founders_given)
		{
			founder_sequence_generator gen;
			output_sequences(gen, args_info);
		}
		else if (args_info.output_haplotypes_given)
		{
			haplotype_sequence_generator gen;
			output_sequences(gen, args_info);
		}
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
