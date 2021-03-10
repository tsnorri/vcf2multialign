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
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "sequence_generator_base.hh"
#include "utility.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;
namespace vgs	= vcf2multialign::variant_graphs;


namespace vcf2multialign {
	
	void direct_matching_sequence_generator::open_output_file(char const *path, output_stream_type &of, lb::writing_open_mode const mode) const
	{
		lb::open_file_for_writing("REF", of, mode);
	}
	
	
	void sequence_generator_base::read_variant_graph(char const *input_graph_path)
	{
		lb::file_istream input_graph_stream;
		lb::open_file_for_reading(input_graph_path, input_graph_stream);
		
		cereal::PortableBinaryInputArchive iarchive(input_graph_stream);
		iarchive(m_graph);
	}
	
	
	void direct_matching_sequence_generator::output_sequences()
	{
		try
		{
			// Setup the progress indicator.
			this->install_progress_indicator();
			
			// Create a string view from the reference.
			std::string_view const reference_sv(m_reference.data(), m_reference.size());
		
			// Output in chunks.
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(m_may_overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			
			auto const stream_count(this->get_stream_count());
			dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << stream_count << " sequences will be written.\n"; });
			std::size_t const chunk_count(std::ceil(1.0 * stream_count / m_chunk_size));
			for (auto const &pair : ranges::view::closed_iota(std::size_t(0), chunk_count) | ranges::view::sliding(2))
			{
				auto const lhs(pair[0]);
				auto const rhs(pair[1]);
				auto const lhsc(m_chunk_size * lhs);
				auto const rhsc(std::min(stream_count, m_chunk_size * rhs));
				dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Processing chunk " << rhs << '/' << chunk_count << "…\n"; });
				output_stream_vector output_files(rhsc - lhsc); // Cannot reuse b.c. lb::file_ostream has a deleted copy constructor.
				
				// Open the output files.
				// The last stream will be REF.
				if (m_output_reference && rhsc == stream_count)
				{
					open_output_file("REF", output_files.back(), mode);
					for (auto &&[i, of] : ranges::view::zip(ranges::view::iota(lhsc, rhsc), output_files) | ranges::view::drop_last(1))
						open_output_file(i, of, mode);
				}
				else
				{
					for (auto &&[i, of] : ranges::view::zip(ranges::view::iota(lhsc, rhsc), output_files))
						open_output_file(i, of, mode);
				}
				
				output_chunk(reference_sv, output_files);
			}
			
			dispatch_async(dispatch_get_main_queue(), ^{
				lb::log_time(std::cerr);
				std::cerr << "Done.\n"; // FIXME: log statistics?
				this->finish_mt();
			});
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			this->log_assertion_failure_and_exit(exc);
		}
		catch (std::exception const &exc)
		{
			this->log_exception_and_exit(exc);
		}
		catch (...)
		{
			this->log_unknown_exception_and_exit();
		}
	}
	
	
	void direct_matching_sequence_generator::output_chunk(std::string_view const &reference_sv, output_stream_vector &output_files)
	{
		// Use a simple queue for keeping track of the aligned positions of the founders.
		// first: node index, second: stream number.
		stream_position_list files_available(output_files.size()), files_waiting;
		
		auto const &ref_positions(m_graph.ref_positions());
		auto const &aln_positions(m_graph.aligned_ref_positions());
		
		{
			typedef vgs::variant_graph::sample_path_vector	sample_paths_type;
			typedef vgs::variant_graph::path_edge_matrix	path_edges_type;
			
			auto handler(this->make_alt_edge_handler(reference_sv, m_graph, output_files));
			auto const &subgraph_start_positions(m_graph.subgraph_start_positions());
			auto const &alt_edge_count_csum(m_graph.alt_edge_count_csum());
			auto const &alt_edge_targets(m_graph.alt_edge_targets());
			auto const &alt_edge_labels(m_graph.alt_edge_labels());
			auto const &sample_paths(m_graph.sample_paths());
			auto const &path_edges(m_graph.path_edges());
			
			// Setup the progress bar.
			progress_indicator_delegate progress_delegate(ref_positions.size() - 2); // FIXME: is the maximum value exact?
			this->progress_indicator().log_with_progress_bar("\t", progress_delegate);
			
			// Setup files_available.
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
							ranges::view::repeat_n(0, (first_subgraph_starts_from_zero ? 0 : 1)),	// Handle the possible padding. (The initial zero could also be included when creating the graph.)
							subgraph_start_positions,
							ranges::view::single(ref_positions.size() - 2)
						) | ranges::view::sliding(2)
					)
				)
			);
			for (auto const &[subgraph_idx, tup] : rsv)
			{
				auto const &[sample_paths, edges_by_path_and_variant, ssp_pair] = tup;
				auto const subgraph_lhs(ssp_pair[0]);
				auto const subgraph_rhs(ssp_pair[1]);
				auto const subgraph_variants(edges_by_path_and_variant.number_of_rows());
				auto const subgraph_paths(edges_by_path_and_variant.number_of_columns());
				std::size_t subgraph_variant_idx(0); // within the subgraph.
				
				handler->set_subgraph_samples_and_edges(sample_paths, edges_by_path_and_variant);
				for (auto const &[node_idx, idx_pair] : ranges::view::enumerate(ranges::view::closed_iota(subgraph_lhs, subgraph_rhs) | ranges::view::sliding(2)))
				{
					auto const lhs(idx_pair[0]);
					auto const rhs(idx_pair[1]);
					auto const alt_lhs(alt_edge_count_csum[lhs]);
					auto const alt_rhs(alt_edge_count_csum[rhs]);
					auto const ref_lhs(ref_positions[1 + lhs]);
					auto const ref_rhs(ref_positions[1 + rhs]);
					auto const aln_lhs(aln_positions[1 + lhs]);
					auto const aln_rhs(aln_positions[1 + rhs]);
					libbio_assert_lte(ref_lhs, ref_rhs);
					libbio_assert_lte(aln_lhs, aln_rhs);
					auto const ref_len(ref_rhs - ref_lhs);
					auto const aln_len(aln_rhs - aln_lhs);
					libbio_assert_lte(ref_len, aln_len);
					
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
						auto const ref_sub(reference_sv.substr(ref_lhs, ref_len));
						for (auto &sp : files_waiting)
						{
							auto &stream(output_files[sp.stream_number]);
							libbio_assert_eq(lhs, sp.node);
							libbio_assert_eq(aln_lhs, stream.tellp());
							stream << ref_sub;
							std::fill_n(std::ostream_iterator <char>(stream), aln_len - ref_len, '-');
							libbio_assert_eq(aln_rhs, stream.tellp());
							sp.node = rhs;
						}
					}
					else
					{
						if (!files_waiting.empty())
						{
							handler->set_position(lhs, rhs, alt_lhs, ref_lhs, ref_rhs, aln_lhs, aln_rhs);
							handler->handle_node_with_alt_edges(files_waiting, subgraph_variant_idx);
						}
						++subgraph_variant_idx;
					}
					
					// Sort the stream numbers by the writing position, i.e. node number.
					files_waiting.sort();
					files_available.merge(files_waiting);
					
					progress_delegate.advance();
				}
			}
			this->end_logging();
		}
		
		// Fill with the reference up to aligned reference length.
		{
			dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Filling with the reference…\n"; });
			
			progress_indicator_delegate progress_delegate(files_available.size());
			this->progress_indicator().log_with_progress_bar("\t", progress_delegate);
			
			auto const ref_end(ref_positions.back());
			for (auto const &sp : files_available)
			{
				// Fill.
				auto const ref_begin(ref_positions[1 + sp.node]);
				auto const ref_sub(reference_sv.substr(ref_begin, ref_end - ref_begin));
				auto &stream(output_files[sp.stream_number]);
				stream << ref_sub;
				stream.flush();
				
				progress_delegate.advance();
			}
			
			this->end_logging();
		}
		
		dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Finishing…\n"; });
	}
	
	
	void alt_edge_handler_base::handle_edge(stream_position &sp, std::size_t const edge_number) const
	{
		auto &stream((*m_output_files)[sp.stream_number]);
		libbio_assert_eq(m_lhs, sp.node);
		libbio_assert_eq(m_aln_lhs, stream.tellp());
		if (0 == edge_number)
		{
			// REF edge.
			auto const aln_len(m_aln_rhs - m_aln_lhs);
			auto const ref_len(m_ref_rhs - m_ref_lhs);
			auto const ref_sub(m_reference_sv->substr(m_ref_lhs, ref_len));
			libbio_assert_lte(ref_len, aln_len);
			stream << ref_sub;
			std::fill_n(std::ostream_iterator <char>(stream), aln_len - ref_len, '-');
			libbio_assert_eq(m_aln_rhs, stream.tellp());
			sp.node = m_rhs;
		}
		else
		{
			// ALT edge.
			auto const &aln_positions(m_graph->aligned_ref_positions());
			auto const &alt_edge_labels(m_graph->alt_edge_labels());
			auto const &alt_edge_targets(m_graph->alt_edge_targets());
			auto const alt_idx(m_alt_lhs + edge_number - 1);
			auto const target_node(alt_edge_targets[alt_idx]);
			auto const &alt(alt_edge_labels[alt_idx]);

			auto const aln_rhs(aln_positions[1 + target_node]);
			libbio_assert_lte(m_aln_lhs, aln_rhs);
			auto const aln_len(aln_rhs - m_aln_lhs);
			libbio_always_assert_lte(alt.size(), aln_len);
				
			stream << alt;
			std::fill_n(std::ostream_iterator <char>(stream), aln_len - alt.size(), '-');
			libbio_assert_eq(aln_rhs, stream.tellp());
			sp.node = target_node;
		}
	}
}
