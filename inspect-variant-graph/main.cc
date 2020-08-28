/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <charconv>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;
namespace vgs	= vcf2multialign::variant_graphs;


namespace {
	
	enum class reading_mode
	{
		NODE_NUMBER,
		REF_POSITION
	};
	
	
	reading_mode g_mode{reading_mode::NODE_NUMBER};
	

	template <typename t_range>
	std::size_t read_next_node_number(t_range const &range)
	{
		auto const length(ranges::size(range));
		std::string buffer;
		
		while (true)
		{
			switch (g_mode)
			{
				case reading_mode::NODE_NUMBER:
				{
					std::cout << "Node number [0, " << length << ") or ‘r’ for REF position? " << std::flush;
					std::cin >> buffer;
					if (std::cin.eof())
						return SIZE_MAX;
					
					if ("r" == buffer)
					{
						g_mode = reading_mode::REF_POSITION;
						continue;
					}

					std::size_t node_number{};
					auto const res(std::from_chars(buffer.data(), buffer.data() + buffer.size(), node_number));
					if (std::errc{} != res.ec)
						continue;
					
					if (! (node_number < length))
						continue;
					
					return node_number;
				}
					
				case reading_mode::REF_POSITION:
				{
					std::cout << "REF position or ‘n’ for node number? " << std::flush;
					std::cin >> buffer;
					if (std::cin.eof())
						return SIZE_MAX;
					
					if ("n" == buffer)
					{
						g_mode = reading_mode::NODE_NUMBER;
						continue;
					}
					
					std::size_t ref_pos{};
					auto const res(std::from_chars(buffer.data(), buffer.data() + buffer.size(), ref_pos));
					if (std::errc{} != res.ec)
						continue;
					
					auto const it(ranges::upper_bound(
						range,
						ref_pos,
						ranges::less{},
						[](auto const &item){
							auto const &[ref_pair, aligned_ref_pair, alt_edge_pair] = item;
							auto const ref_lhs(ref_pair[0]);
							return ref_lhs;
						}
					));
					
					return std::distance(ranges::begin(range), it - 1);
				}
			}
		}
	}
	
	
	void inspect_variant_graph(
		char const *reference_path,
		char const *input_graph_path,
		char const *reference_seq_name
	)
	{
		v2m::vector_type reference;
		vgs::variant_graph graph;
		
		{
			lb::file_istream input_graph_stream;
			
			// Open the files.
			lb::open_file_for_reading(input_graph_path, input_graph_stream);
			
			lb::mmap_handle <char> ref_handle;
			ref_handle.open(reference_path);
			
			// Read the input FASTA.
			v2m::read_single_fasta_seq(ref_handle, reference, reference_seq_name);
			
			// Read the intermediate graph.
			cereal::PortableBinaryInputArchive iarchive(input_graph_stream);
			iarchive(graph);
		}
		
		{
			auto const &ref_positions(graph.ref_positions());
			auto const &aligned_ref_positions(graph.aligned_ref_positions());
			auto const &subgraph_start_positions(graph.subgraph_start_positions());
			auto const &alt_edge_count_csum(graph.alt_edge_count_csum());
			auto const &alt_edge_targets(graph.alt_edge_targets());
			auto const &alt_edge_labels(graph.alt_edge_labels());
			
			auto const rsv(ranges::view::zip(
				ref_positions | ranges::view::tail | ranges::view::sliding(2),
				aligned_ref_positions | ranges::view::tail | ranges::view::sliding(2),
				alt_edge_count_csum | ranges::view::sliding(2)
			));
			
			if (std::cin.eof())
				return;
			while (true)
			{
				auto const idx(read_next_node_number(rsv));
				if (SIZE_MAX == idx)
					break;
				
				auto const &pairs(rsv[idx]);
				auto const &[ref_pair, aligned_ref_pair, alt_edge_pair] = pairs;
				auto const aligned_lhs(aligned_ref_pair[0]);
				auto const aligned_rhs(aligned_ref_pair[1]);
				auto const ref_lhs(ref_pair[0]);
				auto const ref_rhs(ref_pair[1]);
				auto const alt_edge_start(alt_edge_pair[0]);
				auto const alt_edge_limit(alt_edge_pair[1]);
				std::string_view const ref_sub(reference.data() + ref_lhs, ref_rhs - ref_lhs);
				
				std::cout << "REF position:\t\t" << ref_lhs << '\n';
				std::cout << "Next REF position:\t" << ref_rhs << '\n';
				std::cout << "Aligned position:\t" << aligned_lhs << '\n';
				std::cout << "Next aligned position:\t" << aligned_rhs << '\n';
				
				std::cout << "REF (" << ref_sub.size() << "):\t\t";
				if (20 < ref_sub.size())
				{
					std::cout << ref_sub.substr(0, 10);
					std::cout << "…";
					std::cout << ref_sub.substr(ref_sub.size() - 10, 10);
				}
				else
				{
					std::cout << ref_sub;
				}
				std::cout << '\n';

				auto const alt_rsv(ranges::view::slice(ranges::view::zip(alt_edge_targets, alt_edge_labels), alt_edge_start, alt_edge_limit));
				auto const alt_edge_count(ranges::size(alt_rsv));
				std::cout << "ALT edges:\t\t" << alt_edge_count << '\n';
				if (alt_edge_count)
				{
					std::cout << "ALT target          Label\n";
					for (auto const &[target, label] : alt_rsv)
						std::cout << std::left << std::setw(20) << target << label << " (" << label.size() << ")\n";
				}
				std::cout << '\n';
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
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	
	try
	{
		inspect_variant_graph(
			args_info.reference_arg,
			args_info.variants_arg,
			args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr
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
