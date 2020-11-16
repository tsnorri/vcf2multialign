/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */


#include <range/v3/view/subrange.hpp>
#include <vcf2multialign/path_mapping/path_mapper.hh>
#include <vcf2multialign/path_mapping/segment_connector.hh>
#include "founder_sequence_greedy_generator.hh"
#include "utility.hh"

// Substring refers here to a labelled graph path through one subgraph whereas path refers to a graph path through two consecutive subgraphs.

namespace lb	= libbio;
namespace pm	= vcf2multialign::path_mapping;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;
namespace vgs	= vcf2multialign::variant_graphs;


namespace {
	
	class sequence_writer
	{
	public:
		typedef v2m::founder_sequence_greedy_generator::removed_count_map	removed_count_map;
		
	protected:
		std::string_view					m_reference;
		lb::int_vector <1>					m_seen_edges;
		vgs::variant_graph const			*m_graph{};
		removed_count_map					*m_removed_counts{};
		std::size_t							m_tail_length{};
		
	public:
		sequence_writer(
			vgs::variant_graph const &graph,
			v2m::vector_type const &reference,
			removed_count_map &removed_counts,
			std::size_t const tail_length
		):
			m_reference(reference.data(), reference.size()),
			m_graph(&graph),
			m_removed_counts(&removed_counts),
			m_tail_length(tail_length)
		{
		}
		
		void output_ref(std::size_t const node_idx, bool const should_remove_mid, std::ostream &os) const { output_ref(node_idx, 1 + node_idx, should_remove_mid, os); }
		void output_subgraph_path(std::size_t const subgraph_idx, pm::substring_index_type const path_idx, bool const should_remove_mid, std::ostream &os);
		void output_ref_for_subgraph(std::size_t const subgraph_idx, bool const should_remove_mid, std::ostream &os) const;
		void output_gaps_for_subgraph(std::size_t const subgraph_idx, std::ostream &os) const { output_char_for_subgraph(subgraph_idx, '-', os); }
		void output_n_for_subgraph(std::size_t const subgraph_idx, std::ostream &os) const { output_char_for_subgraph(subgraph_idx, 'N', os); }
		
	protected:
		void output_char_for_subgraph(std::size_t const subgraph_idx, char const c, std::ostream &os) const;
		
		// Return true if middle part was removed.
		bool output_ref(pm::substring_index_type const node_idx, pm::substring_index_type const next_node_idx, bool const should_remove_mid, std::ostream &os) const;
		bool output_alt(pm::substring_index_type const node_idx, pm::substring_index_type const next_node_idx, std::size_t const alt_idx, bool const should_remove_mid, std::ostream &os) const;
		bool output_part(std::string_view const &part, std::size_t const gap_count, bool const should_remove_mid, std::ostream &os) const;
		
		// (Ab)use the fact that m_removed_counts is a pointer, which is not modified.
		void mark_middle_part_removed(std::size_t const node_idx, bool const is_ref) const { ++(*m_removed_counts)[node_idx]; }
	};
	
	
	void sequence_writer::output_subgraph_path(
		std::size_t const subgraph_idx,
		pm::substring_index_type const path_idx,
		bool const should_remove_mid,
		std::ostream &os
	)
	{
		auto const &ref_positions(m_graph->ref_positions());						// REF positions (0-based) by node number. We use 1-based indexing in order to make summing easier.
		auto const &aln_positions(m_graph->aligned_ref_positions());				// Aligned REF positions by node number. We use 1-based indexing in order to make summing easier.
		auto const &subgraph_start_positions(m_graph->subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		auto const &alt_edge_targets(m_graph->alt_edge_targets());					// ALT edge target nodes by edge number as a concatenated vector.
		auto const &alt_edge_count_csum(m_graph->alt_edge_count_csum());			// Cumulative sum of ALT edges by node number.
		auto const &all_path_edges(m_graph->path_edges());							// Edge numbers (0 for REF edge, 1 for first ALT edge etc.) by path, variant and subgraph number.
		auto const &subgraph_path_edges(all_path_edges[subgraph_idx]);
		// Get a slice that represents the current path. Indices represent the variants (i.e. nodes that have ALT edges).
		auto const &current_path_edges(subgraph_path_edges.column(path_idx));
		
		// Node indices.
		auto const subgraph_begin(subgraph_start_positions[subgraph_idx]);
		auto const subgraph_end(
			1 + subgraph_idx < subgraph_start_positions.size()
			? subgraph_start_positions[1 + subgraph_idx]
			: ref_positions.size() - 2
		);
		
		libbio_assert_eq(os.tellp(), aln_positions[1 + subgraph_begin]);
		
		std::size_t variant_idx(0);
		auto expected_node_idx(subgraph_begin);
		for (auto node_idx(subgraph_begin); node_idx < subgraph_end; ++node_idx)
		{
			std::fill(m_seen_edges.word_begin(), m_seen_edges.word_end(), 0);
			
			libbio_assert_lt(1 + node_idx, alt_edge_count_csum.size());
			auto const alt_edge_start(alt_edge_count_csum[node_idx]);
			auto const alt_edge_limit(alt_edge_count_csum[1 + node_idx]);
			
			if (alt_edge_start == alt_edge_limit)
			{
				// Non-variant node.
				if (node_idx == expected_node_idx)
				{
					++expected_node_idx;
					if (output_ref(node_idx, expected_node_idx, should_remove_mid, os))
						mark_middle_part_removed(node_idx, true);
				}
			}
			else
			{
				// Variant node.
				if (node_idx == expected_node_idx)
				{
					auto const edge(current_path_edges[variant_idx]);
					if (0 == edge)
					{
						++expected_node_idx;
						if (output_ref(node_idx, expected_node_idx, should_remove_mid, os))
							mark_middle_part_removed(node_idx, true);
					}
					else
					{
						auto const alt_idx(alt_edge_start + edge - 1);
						expected_node_idx = alt_edge_targets[alt_idx];
						libbio_assert_lte_msg(
							expected_node_idx, subgraph_end,
							"Subgraph ", subgraph_idx, ", node ", node_idx, ", alt edge ", edge, ": target node is ",
							alt_edge_targets[alt_idx], " but subgraph end is at node ", subgraph_end
						);
							
						if (! (alt_idx < m_seen_edges.size()))
							m_seen_edges.resize(1 + alt_idx, 0x0);
						
						if (output_alt(node_idx, expected_node_idx, alt_idx, should_remove_mid && !m_seen_edges[alt_idx], os))
							mark_middle_part_removed(node_idx, false);
						else
							m_seen_edges[alt_idx] |= 0x1;
					}
				}
				++variant_idx;
			}
		}
		
		libbio_assert_eq(os.tellp(), aln_positions[1 + subgraph_end]);
	}
	
	
	void sequence_writer::output_ref_for_subgraph(
		std::size_t const subgraph_idx,
		bool const should_remove_mid,
		std::ostream &os
	) const
	{
		auto const &ref_positions(m_graph->ref_positions());
		auto const &subgraph_start_positions(m_graph->subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		
		// Node indices.
		auto const subgraph_begin(subgraph_start_positions[subgraph_idx]);
		auto const subgraph_end(
			1 + subgraph_idx < subgraph_start_positions.size()
			? subgraph_start_positions[1 + subgraph_idx]
			: ref_positions.size() - 2
		);
		
		libbio_assert_eq(os.tellp(), m_graph->aligned_ref_positions()[1 + subgraph_begin]);
		for (auto i(subgraph_begin); i < subgraph_end; ++i)
		{
			if (output_ref(i, 1 + i, should_remove_mid, os))
				mark_middle_part_removed(i, true);
		}
		libbio_assert_eq(os.tellp(), m_graph->aligned_ref_positions()[1 + subgraph_end]);
	}
	
	
	void sequence_writer::output_char_for_subgraph(std::size_t const subgraph_idx, char const cc, std::ostream &os) const
	{
		auto const &subgraph_start_positions(m_graph->subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		auto const &ref_positions(m_graph->ref_positions());
		auto const &aln_positions(m_graph->aligned_ref_positions());
		
		// Node indices.
		auto const subgraph_begin(subgraph_start_positions[subgraph_idx]);
		auto const subgraph_end(
			1 + subgraph_idx < subgraph_start_positions.size()
			? subgraph_start_positions[1 + subgraph_idx]
			: ref_positions.size() - 2
		);
		
		auto const lhs_aln_pos(aln_positions[1 + subgraph_begin]);
		auto const rhs_aln_pos(aln_positions[1 + subgraph_end]);
		libbio_assert_lte(lhs_aln_pos, rhs_aln_pos);
		auto const gap_count(rhs_aln_pos - lhs_aln_pos);
		libbio_assert_eq(os.tellp(), lhs_aln_pos);
		std::fill_n(std::ostream_iterator <char>(os), gap_count, cc);
		libbio_assert_eq(os.tellp(), rhs_aln_pos);
	}
	
	
	bool sequence_writer::output_ref(
		pm::substring_index_type const node_idx,
		pm::substring_index_type const next_node_idx,
		bool const should_remove_mid,
		std::ostream &os
	) const
	{
		auto const &ref_positions(m_graph->ref_positions());
		auto const &aln_positions(m_graph->aligned_ref_positions());
		auto const lhs_ref_pos(ref_positions[1 + node_idx]);
		auto const rhs_ref_pos(ref_positions[1 + next_node_idx]);
		auto const lhs_aln_pos(aln_positions[1 + node_idx]);
		auto const rhs_aln_pos(aln_positions[1 + next_node_idx]);
		libbio_assert_lte(lhs_ref_pos, rhs_ref_pos);
		libbio_assert_lte(lhs_aln_pos, rhs_aln_pos);
		auto const ref_len(rhs_ref_pos - lhs_ref_pos);
		auto const aln_len(rhs_aln_pos - lhs_aln_pos);
		libbio_assert_lte(ref_len, aln_len);
		auto const gap_count(aln_len - ref_len);
		auto const ref_sub(m_reference.substr(lhs_ref_pos, ref_len));
		
		return output_part(ref_sub, gap_count, should_remove_mid, os);
	}
	
	
	bool sequence_writer::output_alt(
		pm::substring_index_type const node_idx,
		pm::substring_index_type const next_node_idx,
		std::size_t const alt_idx,
		bool const should_remove_mid,
		std::ostream &os
	) const
	{
		auto const &aln_positions(m_graph->aligned_ref_positions());
		auto const &alt_edge_labels(m_graph->alt_edge_labels());
		auto const lhs_aln_pos(aln_positions[1 + node_idx]);
		auto const rhs_aln_pos(aln_positions[1 + next_node_idx]);
		libbio_assert_lte(lhs_aln_pos, rhs_aln_pos);
		auto const alt_str(alt_edge_labels[alt_idx]);
		auto const aln_len(rhs_aln_pos - lhs_aln_pos);
		libbio_assert_lte(alt_str.size(), aln_len);
		auto const gap_count(aln_len - alt_str.size());
		
		bool const retval(output_part(alt_str, gap_count, should_remove_mid, os));
		std::fill_n(std::ostream_iterator <char>(os), gap_count, '-');
		
		return retval;
	}
	
	
	bool sequence_writer::output_part(
		std::string_view const &part,
		std::size_t const gap_count,
		bool const should_remove_mid,
		std::ostream &os
	) const
	{
		if (should_remove_mid && 2 * m_tail_length <= part.size())
		{
			auto const head(part.substr(0, m_tail_length));
			auto const tail(part.substr(part.size() - m_tail_length));
			auto const mid_len(part.size() - 2 * m_tail_length);
			os << head;
			std::fill_n(std::ostream_iterator <char>(os), mid_len, '-');
			os << tail;
			return true;
		}
		else
		{
			os << part;
			std::fill_n(std::ostream_iterator <char>(os), gap_count, '-');
			return false;
		}
	}
}


namespace vcf2multialign {
	
	void founder_sequence_greedy_generator::process_graph_and_output(
		output_stream_vector &output_files,
		removed_count_map &removed_counts,
		progress_indicator_delegate &progress_delegate
	) const
	{
		// Greedy matching is done as follows:
		// 1. Initially, occurring substring numbers in the first segment are assigned to each slot (i.e. output stream).
		//    The remaining slots are marked “free”. (Done with slots_available_lhs.)
		// 2. For each segment pair, the right hand side slots are treated similarly. (Done with slots_available_rhs.)
		// 3. The sorted (by count) path list is then traversed. If a pair of substring numbers is available, they are chosen.
		//    If either lhs or rhs is available and the other side has a free slot, that slot is assigned the corresponding substring index.
		//    If neither lhs or rhs is available but both sides have a free slot, the slots are assigned the corresponding substring indices.
		// 4. Finally, the edges are drawn and the remaining substrings are connected arbitrarily.
		// 5. If a slot is available after all edges have been connected, it is returned to the free list and filled with gap characters.
		
		auto const &sample_paths(m_graph.sample_paths());							// Sample path numbers by sample and subgraph number, vector of vectors.
		auto const &all_path_edges(m_graph.path_edges());
		auto const &subgraph_start_positions(m_graph.subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		
		// Before beginning, check whether the first subgraph starts from zero.
		if (subgraph_start_positions.empty())
			return;
		
		sequence_writer sw(m_graph, m_reference, removed_counts, m_tail_length);
		
		bool const first_subgraph_starts_from_zero(0 == subgraph_start_positions.front());
		if (!first_subgraph_starts_from_zero)
		{
			// Output REF until the starting position.
			auto const first_subgraph_start_pos(subgraph_start_positions.front());
			for (auto &output_file : output_files)
			{
				for (std::size_t i(0); i < first_subgraph_start_pos; ++i)
					sw.output_ref(i, m_should_remove_mid, output_file);
			}
		}
		
		// Count the distinct pairs of paths.
		// A vector is used b.c. we would like to sort by different keys.
		
		pm::path_item_vector path_counts;					// Existing edges in the graph, to be sorted by occurrence count.
		pm::edge_vector edges;								// Edges to be drawn in the founders.
		pm::substring_index_vector substrings_added_to_lhs;	// Substrings added to the lhs segment in the current iteration.
		pm::segment_connector sc(m_founder_count);
		pm::path_mapper pm(m_founder_count);
		
		// Mark the indices in the first segment as assigned for greedy matching.
		std::size_t max_lhs_substring_idx(all_path_edges.front().number_of_columns() - 1);
		libbio_always_assert_lt_msg(max_lhs_substring_idx, m_founder_count, "Given founder count (", m_founder_count, ") is less than the number of distinct substrings in subgraph 1 (", 1 + max_lhs_substring_idx, ").");
		sc.setup(1 + max_lhs_substring_idx);
		pm.setup(1 + max_lhs_substring_idx);
		
		// Handle the segment pairs.
		for (auto const &[lhs_subgraph_idx, pair] : rsv::zip(rsv::ints(0), sample_paths | rsv::sliding(2)))
		{
			try
			{
				// Initialize for the current iteration.
				path_counts.clear();
				edges.clear();
				substrings_added_to_lhs.clear();
				
				// Destructure the pair.
				auto const &lhs_substring_numbers(pair[0]);
				auto const &rhs_substring_numbers(pair[1]);
				libbio_always_assert_eq(lhs_substring_numbers.size(), rhs_substring_numbers.size()); // Verify that the input is valid.
				
				// Iterate over the pairs of paths and count.
				pm::substring_index_type max_rhs_substring_idx(0);
				for (auto const &[lhs_substring_idx, rhs_substring_idx] : rsv::zip(lhs_substring_numbers, rhs_substring_numbers))
				{
					// Try to find an existing path pair.
					auto const res(std::equal_range(path_counts.begin(), path_counts.end(), pm::path_item(lhs_substring_idx, rhs_substring_idx)));
					if (res.first == res.second)
					{
						// Not found (since there were no not-less-than elements). Append to the end and rotate the greater-than range.
						auto const first_greater_than_idx(std::distance(path_counts.begin(), res.second)); // Index of the first greater-than element.
						path_counts.emplace_back(lhs_substring_idx, rhs_substring_idx, 1); // Invalidates iterators.
						
						libbio_assert_lt(first_greater_than_idx, path_counts.size());
						auto const begin(path_counts.begin() + first_greater_than_idx);
						auto const end(path_counts.end());
						std::rotate(begin, end - 1, end); // There is at least one element b.c. we just emplace_back’d one.
					}
					else
					{
						// Found.
						libbio_assert_eq(1, std::distance(res.first, res.second)); // Check that there is in fact only one matching item.
						res.first->count += 1;
					}
					
					max_rhs_substring_idx = lb::max_ct(max_rhs_substring_idx, rhs_substring_idx);
				}
				
				libbio_assert(std::is_sorted(path_counts.begin(), path_counts.end()));
				libbio_always_assert_lt_msg(max_rhs_substring_idx, m_founder_count, "Given founder count (", m_founder_count, ") is less than the number of distinct substrings in subgraph ", lhs_subgraph_idx, " (", 1 + max_rhs_substring_idx, ").");
				
				// Sort by count and add edges in descending count order.
				std::sort(path_counts.begin(), path_counts.end(), [](auto const &lhs, auto const &rhs) -> bool {
					return lhs.count < rhs.count;
				});
				sc.make_edges(path_counts, 1 + max_rhs_substring_idx, edges, substrings_added_to_lhs);
				
				// Associate the edges with the founders i.e. output streams.
				// assign_edges_to_founders requires the edges to be sorted by lhs_idx.
				std::sort(edges.begin(), edges.end());
				pm.add_substrings(substrings_added_to_lhs);
				pm.assign_edges_to_founders(edges);
				pm.update_string_indices();
				
				// Output the founders.
				for (auto const &[founder_idx, substring_idx] : rsv::enumerate(pm.string_indices_by_founder()))
				{
					libbio_assert_msg(!m_output_reference || founder_idx != output_files.size() - 1, "Expected founder sequence not to be output to the reference file");
					auto &os(output_files[founder_idx]);
					if (pm::UNASSIGNED_INDEX == substring_idx)
					{
						if (m_fill_unassigned_with_ref)
							sw.output_ref_for_subgraph(lhs_subgraph_idx, m_should_remove_mid, os);
						else
							sw.output_n_for_subgraph(lhs_subgraph_idx, os);
					}
					else
					{
						sw.output_subgraph_path(lhs_subgraph_idx, substring_idx, m_should_remove_mid, os);
					}
				}
				
				if (m_output_reference)
				{
					auto &os(output_files.back());
					sw.output_ref_for_subgraph(lhs_subgraph_idx, false, os);
				}
				
				pm.end_subgraph();
				progress_delegate.advance();
				max_lhs_substring_idx = max_rhs_substring_idx;
			}
			catch (lb::assertion_failure_exception const &)
			{
				std::cerr << "Assertion failure in subgraph " << lhs_subgraph_idx << ":\n";
				std::cerr << "substrings_added_to_lhs:";
				for (auto const substring_idx : substrings_added_to_lhs)
					std::cerr << ' ' << substring_idx;
				std::cerr << '\n';
				throw;
			}
		}

		// Handle the last subgraph and output the founders.
		{
			pm.update_string_indices();
			auto const last_subgraph_idx(m_graph.subgraph_count() - 1);
			for (auto const &[founder_idx, substring_idx] : rsv::enumerate(pm.string_indices_by_founder()))
			{
				auto &os(output_files[founder_idx]);
				if (pm::UNASSIGNED_INDEX == substring_idx)
				{
					if (m_fill_unassigned_with_ref)
						sw.output_ref_for_subgraph(last_subgraph_idx, m_should_remove_mid, os);
					else
						sw.output_n_for_subgraph(last_subgraph_idx, os);
				}
				else
				{
					sw.output_subgraph_path(last_subgraph_idx, substring_idx, m_should_remove_mid, os);
				}
				os << std::flush;
			}
			
			if (m_output_reference)
			{
				auto &os(output_files.back());
				sw.output_ref_for_subgraph(last_subgraph_idx, false, os);
				os << std::flush;
			}
		}
		progress_delegate.advance();
	}
	
	
	void founder_sequence_greedy_generator::output_sequences()
	{
		try
		{
			// Setup the progress indicator.
			this->install_progress_indicator();
			progress_indicator_delegate progress_delegate(m_graph.subgraph_count());
			this->progress_indicator().log_with_progress_bar("\t", progress_delegate);
			
			// Create a string view from the reference.
			std::string_view const reference_sv(m_reference.data(), m_reference.size());
			
			auto const output_count(m_founder_count + m_output_reference);
			dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << m_founder_count << " sequences will be written.\n"; });
			
			// Don’t use chunks for now b.c. that would complicate things too much and not writing everything simultaneously
			// is more useful with predicted sequences (not founders).
			
			output_stream_vector output_files(output_count);
			
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(m_may_overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			
			for (std::size_t i(0); i < m_founder_count; ++i)
				open_founder_output_file(i, output_files[i], mode);
			
			if (m_output_reference)
				lb::open_file_for_writing("REF", output_files.back(), mode);
			
			removed_count_map removed_counts;
			
			// Generate the sequences.
			process_graph_and_output(output_files, removed_counts, progress_delegate);
			
			dispatch_async(dispatch_get_main_queue(), ^{
				std::cerr << '\n';
				lb::log_time(std::cerr);
				std::cerr << "Done.\n";
				
				{
					std::cout << "# Nodes having edge labels the middle part of which was removed:\n";
					for (auto const &[node_idx, removed_count] : removed_counts)
						std::cout << node_idx << '\t' << removed_count << '\n';
				}
				
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
}
