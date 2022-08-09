/*
 * Copyright (c) 2019–2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <unistd.h>
#include <vcf2multialign/utility/dispatch_exit_guard.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include "cmdline.h"
#include "founder_sequence_greedy_generator.hh"
#include "sequence_generator_base.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;
namespace vgs	= vcf2multialign::variant_graphs;


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
				auto const path_number((*m_sample_paths)[sp.stream_number]);
				auto const edge_number((*m_edges_by_path_and_variant)(subgraph_variant_idx, path_number));
				this->handle_edge(sp, edge_number);
			}
		}
	};
	
	
	class founder_sequence_generator final : public v2m::direct_matching_sequence_generator
	{
	public:
		typedef v2m::direct_matching_sequence_generator::output_stream_vector	output_stream_vector;
		
		using v2m::direct_matching_sequence_generator::direct_matching_sequence_generator;
		
		std::unique_ptr <v2m::alt_edge_handler_base> make_alt_edge_handler(
			std::string_view const &reference_sv,
			vgs::variant_graph const &graph,
			output_stream_vector &output_files
		) const override
		{
			return std::unique_ptr <v2m::alt_edge_handler_base>(new founder_alt_edge_handler(reference_sv, graph, output_files));
		}
		
		std::size_t get_stream_count() const override
		{
			auto const max_paths(m_graph.max_paths_in_subgraph());
			return (m_output_reference ? 1 : 0) + max_paths;
		}
		
		std::string output_path(std::size_t const file_idx) const override
		{
			return std::to_string(file_idx);
		}
	};
	
	
	class sample_sequence_generator final : public v2m::direct_matching_sequence_generator
	{
	protected:
		// FIXME: We currently assume that all samples (and records) have the same ploidy.
		std::size_t m_sample_count{};	// Not taking REF into account.
		std::size_t m_ploidy{};
		
	public:
		typedef v2m::direct_matching_sequence_generator::output_stream_vector	output_stream_vector;
		
		using v2m::direct_matching_sequence_generator::direct_matching_sequence_generator;
		
		void prepare_variant_graph() override
		{
			auto const &sample_paths_by_subgraph(m_graph.sample_paths());
			m_sample_count = sample_paths_by_subgraph.empty() ? 0 : (sample_paths_by_subgraph.front().size());
			m_ploidy = m_sample_count / m_graph.sample_names().size();
		}
		
		std::unique_ptr <v2m::alt_edge_handler_base> make_alt_edge_handler(
			std::string_view const &reference_sv,
			vgs::variant_graph const &graph,
			output_stream_vector &output_files
		) const override
		{
			return std::unique_ptr <v2m::alt_edge_handler_base>(new haplotype_alt_edge_handler(reference_sv, graph, output_files));
		}
		
		std::size_t get_stream_count() const override
		{
			return (m_output_reference ? 1 : 0) + m_sample_count;
		}
		
		std::string output_path(std::size_t const file_idx) const override
		{
			auto const name_idx(file_idx / m_ploidy);
			auto const chr_idx(file_idx % m_ploidy + 1);
			auto const &sample_name(m_graph.sample_names()[name_idx]);
			return boost::str(boost::format("%s-%d") % sample_name % chr_idx);
		}
	};
	
	
	template <typename t_generator>
	std::unique_ptr <v2m::dispatch_exit_guard_helper <t_generator>> instantiate_generator(gengetopt_args_info const &args_info)
	{
		return std::make_unique <v2m::dispatch_exit_guard_helper <t_generator>>(
			args_info.pipe_arg,
			args_info.chunk_size_arg,
			(args_info.omit_reference_output_flag ? false : true),
			args_info.overwrite_flag
		);
	}
	
	
	void prepare(v2m::sequence_generator_base &gen, gengetopt_args_info const &args_info)
	{
		gen.read_reference(
			args_info.reference_arg,
			(args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr)
		);
		
		lb::log_time(std::cerr);
		std::cerr << "Reading the variant graph…\n";
		gen.read_variant_graph(args_info.variants_arg);
		gen.prepare_variant_graph();
	}
	
	
	template <typename t_generator>
	void output_sequences(std::unique_ptr <v2m::dispatch_exit_guard_helper <t_generator>> &&gen_ptr)
	{
		// Run in background in order to be able to update a progress bar.
		lb::log_time(std::cerr);
		std::cerr << "Outputting the sequences.\n";
		lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
			[
				gen_ptr = std::move(gen_ptr)
			](){
				gen_ptr->value.run(true);
			}
		);
	}
}


// Try to prevent inlining b.c. dispatch_main() does not work well with a try block in the same function.
// FIXME: some compiler attribute would be better b.c. I don’t think extern actually prevents the call in main() from being inlined.
extern void process(gengetopt_args_info &args_info)
{
	try
	{
		if (args_info.output_founders_given)
		{
			auto gen_ptr(instantiate_generator <founder_sequence_generator>(args_info));
			prepare(gen_ptr->value, args_info);
			output_sequences(std::move(gen_ptr));
		}
		else if (args_info.output_founders_greedy_given)
		{
			typedef v2m::dispatch_exit_guard_helper <v2m::founder_sequence_greedy_generator> wrapped_generator_type;
			auto gen_ptr(std::make_unique <wrapped_generator_type>(
				args_info.pipe_arg,
				args_info.founder_count_arg,
				args_info.tail_length_arg,
				(args_info.omit_reference_output_flag ? false : true),
				args_info.fill_unassigned_with_ref_flag,
				args_info.remove_mid_from_duplicates_flag,
				args_info.overwrite_flag
			));
				
			auto &generator(gen_ptr->value);
			prepare(generator, args_info);
			
			{
				auto const max_paths_in_subgraph(generator.variant_graph().max_paths_in_subgraph());
				if (args_info.founder_count_arg < max_paths_in_subgraph)
				{
					std::cerr
						<< "NOTE: There are at most "
						<< max_paths_in_subgraph
						<< " paths in a subgraph, which is more than the specified number of founder sequences to be generated, "
						<< args_info.founder_count_arg
						<< ".\n";
				}
			}
			
			output_sequences(std::move(gen_ptr));
		}
		else if (args_info.output_samples_given)
		{
			auto gen_ptr(instantiate_generator <sample_sequence_generator>(args_info));
			prepare(gen_ptr->value, args_info);
			output_sequences(std::move(gen_ptr));
		}
	}
	catch (lb::assertion_failure_exception const &exc)
	{
		v2m::log_assertion_failure_exception(exc);
		throw exc;
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
	
	if (args_info.output_founders_greedy_given && args_info.founder_count_arg < 1)
	{
		std::cerr << "Founder count must be positive." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	
	if (args_info.tail_length_arg < 0)
	{
		std::cerr << "Tail length must be non-negative." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	
	if (args_info.founder_count_given && args_info.chunk_size_given)
	{
		std::cerr << "WARNING: Ignoring chunk size when using greedy matching.\n";
	}
	
	{
		auto const mode_count(args_info.output_samples_given + args_info.output_founders_given + args_info.output_founders_greedy_given);
		if (0 == mode_count)
		{
			std::cerr << "No mode given." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	
	// Since we are going to fork, we would like to log the error conditions of the child processes and exit.
	lb::install_dispatch_sigchld_handler(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
		bool did_report_error(false);
		while (true)
		{
			int status(0);
			auto const pid(::waitpid(-1, &status, WNOHANG | WUNTRACED));
			if (pid <= 0)
				break;
			
			if (WIFEXITED(status))
			{
				auto const exit_status(WEXITSTATUS(status));
				if (0 != exit_status)
				{
					did_report_error = true;
					std::cerr << "ERROR: Child process " << pid << " exited with status " << exit_status;
					switch (exit_status)
					{
						case 127:
							std::cerr << " (command not found)";
							break;
							
						case 126:
							std::cerr << " (command invoked cannot execute)";
							break;
						
						case 69: // EX_UNAVAILABLE
							std::cerr << " (service unavailable)";
							break;
						
						case 71: // EX_OSERR
							std::cerr << " (unknown error from execvp())";
							break;
						
						case 74: // EX_IOERR
							std::cerr << " (an I/O error occurred)";
							break;
						
						default:
							break;
					}
					std::cerr << '.' << std::endl;
				}
			}
			else if (WIFSIGNALED(status))
			{
				did_report_error = true;
				auto const signal_number(WTERMSIG(status));
				std::cerr << "ERROR: Child process " << pid << " received signal " << signal_number << '.' << std::endl;
			}
		}
		
		if (did_report_error)
			std::exit(EXIT_FAILURE);
	});

	process(args_info);

	// Needs to be outside the try block in process().
	// Also it seemed that the compiler attempted to optimize the block somehow when
	// it was located in the same function as the call to dispatch_main().
	dispatch_main();

	// Not reached.
	return EXIT_SUCCESS;
}
