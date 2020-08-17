/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/file_handling.hh>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>
#include <set>
#include <vcf2multialign/graph/variant_graph.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	
	void output_substring_counts(v2m::variant_graph const &graph)
	{
		auto const &ref_positions(graph.ref_positions()); // 1-based.
		auto const &aln_positions(graph.aligned_ref_positions()); // 1-based.
		auto const &subgraph_start_positions(graph.subgraph_start_positions());
		auto const &path_edges(graph.path_edges());
		libbio_always_assert_eq(subgraph_start_positions.size(), path_edges.size());
		
		std::size_t max_substrings(0), min_substrings(SIZE_MAX);
		
		std::cout << "SUBGRAPH\tSTART_POSITION\tALIGNED_START_POSITION\tPATHS\n";
		auto const view(rsv::zip(rsv::ints(0), subgraph_start_positions, path_edges));
		for (auto const &[subgraph_idx, start_pos_idx, path_edges] : view)
		{
			auto const ref_pos(ref_positions[1 + start_pos_idx]);
			auto const aln_pos(aln_positions[1 + start_pos_idx]);
			auto const substring_count(path_edges.number_of_columns());
			
			std::cout << subgraph_idx << '\t' << ref_pos << '\t' << aln_pos << '\t' << substring_count << '\n';
			
			max_substrings = std::max(max_substrings, substring_count);
			min_substrings = std::min(min_substrings, substring_count);
		}
		
		std::cerr << "Max substrings: " << max_substrings << '\n';
		std::cerr << "Min substrings: " << min_substrings << '\n';
	}
	
	
	void output_substring_combination_counts(v2m::variant_graph const &graph)
	{
		auto const &sample_paths(graph.sample_paths()); // Sample path numbers by sample and subgraph number, vector of matrices.
		
		typedef std::uint32_t substring_index_type;
		std::set <std::pair <substring_index_type, substring_index_type>> substring_pairs; // Using an std::set is the simplest way to get the count of unique pairs.
		
		std::size_t max_combinations(0), min_combinations(SIZE_MAX);
		
		std::cout << "LHS_SUBGRAPH\tSUBSTRING_COMBINATIONS\n";
		for (auto const &[first_subgraph_idx, pair] : rsv::zip(rsv::ints(0), sample_paths | rsv::sliding(2)))
		{
			auto const &lhs_substring_numbers(pair[0]);
			auto const &rhs_substring_numbers(pair[1]);
			libbio_always_assert_eq(lhs_substring_numbers.size(), rhs_substring_numbers.size());
			
			for (auto const &[lhs_substring_idx, rhs_substring_idx] : rsv::zip(lhs_substring_numbers, rhs_substring_numbers))
			{
				libbio_always_assert_lte(lhs_substring_idx, std::numeric_limits <substring_index_type>::max());
				libbio_always_assert_lte(rhs_substring_idx, std::numeric_limits <substring_index_type>::max());
				substring_pairs.emplace(lhs_substring_idx, rhs_substring_idx);
			}
			
			auto const combination_count(substring_pairs.size());
			
			max_combinations = std::max(max_combinations, combination_count);
			min_combinations = std::min(min_combinations, combination_count);
			
			std::cout << first_subgraph_idx << '\t' << combination_count << '\n';
		}
		
		std::cerr << "Max combinations: " << max_combinations << '\n';
		std::cerr << "Min combinations: " << min_combinations << '\n';
	}


	void output_alt_label_counts(v2m::variant_graph const &graph)
	{
		auto const &alt_labels(graph.alt_edge_labels());
		std::map <std::string, std::size_t> counts;

		for (auto const &label : alt_labels)
			++counts[label];

		for (auto const &[label, count] : counts)
			std::cout << label << '\t' << count << '\n';
	}


	void output_subgraph_lengths(v2m::variant_graph const &graph)
	{
		auto const &subgraph_start_positions(graph.subgraph_start_positions());
		auto const &ref_positions(graph.ref_positions());
		auto const &aligned_ref_positions(graph.aligned_ref_positions());
		auto const view(
			rsv::concat(
				subgraph_start_positions,
				rsv::single(aligned_ref_positions.size() - 2)
			)
		);
		std::cout << "SUBGRAPH_INDEX\tSUBGRAPH_ALIGNED_LENGTH\tSUBGRAPH_REF_LENGTH\tBRIDGE_LENGTH\n";
		for (auto const &[subgraph_idx, pair] : rsv::enumerate(view | rsv::sliding(2)))
		{
			auto const sn(pair[0]);
			auto const en(pair[1]);
			auto const rsp(ref_positions[sn + 1]);
			auto const rep(ref_positions[en + 1]);
			auto const asp(aligned_ref_positions[sn + 1]);
			auto const aep(aligned_ref_positions[en + 1]);
			auto const prev_node_ref_pos(ref_positions[sn]);
			auto const prev_node_aln_pos(aligned_ref_positions[sn]);
			auto const bridge_length(asp - prev_node_aln_pos);
			libbio_always_assert_eq(bridge_length, rsp - prev_node_ref_pos);
			std::cout << subgraph_idx << '\t' << (aep - asp) << '\t' << (rep - rsp) << '\t' << bridge_length << '\n';
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
	
	v2m::variant_graph graph;
	
	{
		lb::file_istream input_graph_stream;
		
		// Open the graph file and read.
		lb::open_file_for_reading(args_info.variants_arg, input_graph_stream);
		cereal::PortableBinaryInputArchive iarchive(input_graph_stream);
		iarchive(graph);
	}
	
	if (args_info.substring_counts_given)
		output_substring_counts(graph);
	else if (args_info.substring_combination_counts_given)
		output_substring_combination_counts(graph);
	else if (args_info.alt_label_counts_given)
		output_alt_label_counts(graph);
	else if (args_info.subgraph_lengths_given)
		output_subgraph_lengths(graph);
	else
	{
		std::cerr << "No mode given.\n";
		std::exit(EXIT_FAILURE);
	}
	
	// Not reached.
	return EXIT_SUCCESS;
}
