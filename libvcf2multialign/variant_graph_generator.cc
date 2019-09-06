/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/algorithm.hh>
#include <libbio/bits.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/transform.hpp>
#include <vcf2multialign/graph/variant_graph_generator.hh>
#include <vcf2multialign/utility/can_handle_variant_alts.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void variant_graph_generator::update_sample_names()
	{
		// XXX I don’t remember why vcf_reader outputs the sample names as a map with names as keys and positions as indices but it should be fixed unless a good reason is found not to.
		auto const &sample_name_map(m_reader->sample_names());
		m_sample_names.resize(sample_name_map.size());
		for (auto const & [sample_name, idx1] : sample_name_map)
			m_sample_names[idx1 - 1] = sample_name; // Copy.
	}
	
	
	void variant_graph_generator::generate_graph()
	{
		update_sample_names();
		m_graph.clear();
		
		m_processed_count.store(0, std::memory_order_relaxed);
		m_reader->reset();
		m_reader->set_parsed_fields(lb::vcf_field::ALL);
		
		auto const &cut_positions(*m_cut_position_list);
		auto const cut_pos_rng(cut_positions.positions | ranges::view::tail);
		auto cut_pos_it(ranges::begin(cut_pos_rng));
		auto const cut_pos_end(ranges::end(cut_pos_rng));
		auto line_no_it(ranges::begin(cut_positions.handled_line_numbers));
		auto const line_no_end(ranges::end(cut_positions.handled_line_numbers));
		std::size_t overlap_end_pos(0);
		std::size_t prev_overlap_end_pos(0);
		libbio_always_assert_neq(cut_pos_it, cut_pos_end);
		
		// Process the variants.
		bool should_continue(false);
		do
		{
			m_reader->fill_buffer();
			bool should_continue_line_number_loop(false);
			should_continue = m_reader->parse(
				[
					this,
					&cut_pos_it,
					cut_pos_end,
					&line_no_it,
					line_no_end,
					&overlap_end_pos,
					&prev_overlap_end_pos
				](lb::transient_variant const &var) -> bool
				{
					auto const lineno(var.lineno());
					libbio_always_assert_lte(lineno, *line_no_it);
					if (lineno == *line_no_it)
					{
						auto const var_pos(var.zero_based_pos());
						auto const var_end(lb::variant_end_pos(var, *m_end_field));
						libbio_always_assert_lte(var_pos, var_end);
						libbio_always_assert_lte(var_pos, *cut_pos_it);
						
						if (var_pos == *cut_pos_it)
						{
							// Handle the subgraph.
							process_subgraph(prev_overlap_end_pos);
							prev_overlap_end_pos = overlap_end_pos;
							
							++cut_pos_it;
							libbio_always_assert_neq(cut_pos_it, cut_pos_end);
						}
						
						// Add to the stack.
						overlap_end_pos = std::max(overlap_end_pos, var_end);
						m_subgraph_variants.emplace_back(var);
						
						++m_processed_count;
						++line_no_it;
						
						auto const should_continue(line_no_it != line_no_end);
						return should_continue;
					}
					return true;
				}
			);
			
		} while (should_continue);
		
		// Process the remaining variants.
		process_subgraph(prev_overlap_end_pos);
		
		// Copy the sample names and finalize.
		m_graph.sample_names() = m_sample_names;
		finalize_graph();
	}
	
	
	// Retrieve the accumulated group of variants and pass them to the worker thread for processing.
	void variant_graph_generator::process_subgraph(std::size_t const prev_overlap_end_pos)
	{
		// Fast path: empty stack.
		if (m_subgraph_variants.empty())
			return;
		
		// Slow path: possibly overlapping variants.
		m_sample_sorter.prepare_for_next_subgraph();
		
		// Sort by variant and ALT.
		for (auto const &var : m_subgraph_variants)
		{
			for (std::size_t i(0), count(var.alts().size()); i < count; ++i)
				m_sample_sorter.sort_by_variant_and_alt(var, 1 + i);
		}
		
		auto const sample_count(m_sample_indexer.total_samples());
		auto const path_count(m_sample_sorter.path_count());
		m_delegate->variant_graph_generator_will_handle_subgraph(m_subgraph_variants.front(), m_subgraph_variants.size(), path_count);
		
		// Begin a new subgraph; place the beginning between the last variant of the previous subgraph and the first variant of the next subgraph.
		std::size_t const subgraph_start_pos(prev_overlap_end_pos + std::ceil((m_subgraph_variants.front().zero_based_pos() - prev_overlap_end_pos) / 2.0));
		auto const [subgraph_start_node_idx, subgraph_start_alt_edge_start_idx, did_create_node_for_subgraph_start] = m_graph.add_main_node(subgraph_start_pos, 0);
		auto const subgraph_idx(m_graph.add_subgraph(subgraph_start_node_idx, sample_count, m_subgraph_variants.size(), path_count));
		if (did_create_node_for_subgraph_start)
			calculate_aligned_ref_pos_for_new_node(subgraph_start_node_idx);
		
		// Iterate over the variants, collect almost everything except for the ALT indices for the paths.
		// The following loop should create exactly one main node with non-zero handled_alt_count per a group of variants with the same POS.
		libbio_always_assert(m_overlap_stack.empty());
		for (auto const &[var_idx, var] : ranges::view::enumerate(m_subgraph_variants))
		{
			auto const ref_pos(var.zero_based_pos());
			
			// Calculate new indices for the handled ALTs.
			auto const total_alt_count(var.alts().size());
			auto const handled_alt_count(std::count_if(var.alts().begin(), var.alts().end(), [](auto const &alt){ return can_handle_variant_alt(alt); }));
			
			// First check if there are ALT edges that need to be handled.
			while (!m_overlap_stack.empty())
			{
				auto const &top_entry(m_overlap_stack.top());
				auto const &prev_var(*top_entry.variant);
				auto const prev_end_pos(lb::variant_end_pos(prev_var, *m_end_field));
				if (ref_pos < prev_end_pos)
				{
					auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(ref_pos, handled_alt_count);
					if (did_create)
						calculate_aligned_ref_pos_for_new_node(node_idx);
					assign_alt_edge_labels_and_queue(var, node_idx, alt_edge_start_idx);
					goto did_handle_current_variant;
				}
				else if (ref_pos == prev_end_pos)
				{
					auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(prev_end_pos, handled_alt_count);
					if (did_create)
						calculate_aligned_ref_pos_for_new_node(node_idx, top_entry.max_alt_edge_aligned_dst_pos);
					m_graph.connect_alt_edges(top_entry.node_number, node_idx);
					assign_alt_edge_labels_and_queue(var, node_idx, alt_edge_start_idx);
					m_overlap_stack.pop();
					goto did_handle_current_variant;
				}
				else // ref_pos > prev_end_pos
				{
					auto const [prev_node_idx, prev_alt_edge_start_idx, did_create] = m_graph.add_main_node(prev_end_pos, 0);
					if (did_create)
						calculate_aligned_ref_pos_for_new_node(prev_node_idx, top_entry.max_alt_edge_aligned_dst_pos);
					m_graph.connect_alt_edges(top_entry.node_number, prev_node_idx);
					m_overlap_stack.pop();
				}
			}
			
			// Handle the current variant if not already done.
			{
				auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(ref_pos, handled_alt_count);
				if (did_create)
					calculate_aligned_ref_pos_for_new_node(node_idx);
				assign_alt_edge_labels_and_queue(var, node_idx, alt_edge_start_idx);
			}
			
		did_handle_current_variant:
			;
		}
		
		// Handle the remaining variants.
		while (!m_overlap_stack.empty())
		{
			auto const &top_entry(m_overlap_stack.top());
			auto const &var(*top_entry.variant);
			auto const end_pos(lb::variant_end_pos(var, *m_end_field));
			auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(end_pos, 0);
			if (did_create)
				calculate_aligned_ref_pos_for_new_node(node_idx, top_entry.max_alt_edge_aligned_dst_pos);
			m_graph.connect_alt_edges(top_entry.node_number, node_idx);
			m_overlap_stack.pop();
		}
		
		// Path number for each sample.
		auto const &paths_by_sample(m_sample_sorter.paths_by_sample());
		auto &paths_by_subgraph(m_graph.sample_paths());
		auto &current_subgraph_paths(paths_by_subgraph[subgraph_idx]);
		for (auto &&[src, dst] : ranges::view::zip(paths_by_sample, current_subgraph_paths))
			dst |= src;
		
		// Get a representative sample for each path and store its fixed ALT indices over the variants with the corresponding path.
		path_sorted_variant psv;
		psv.set_paths_by_sample(m_sample_sorter.paths_by_sample());
		psv.reserve_memory_for_representatives(path_count);
		psv.determine_representatives_for_each_sample();
		auto &dst_path_edges(m_graph.path_edges()[subgraph_idx]);
		std::vector <std::size_t> unhandled_alt_csum;
		for (auto const &[var_idx, var] : ranges::view::enumerate(m_subgraph_variants))
		{
			auto const &alts(var.alts());
			unhandled_alt_csum.clear();
			unhandled_alt_csum.resize(1 + alts.size());
			for (auto const &[i, alt] : ranges::view::enumerate(alts))
				unhandled_alt_csum[1 + i] = unhandled_alt_csum[i] + (can_handle_variant_alt(alt) ? 0 : 1);
			
			auto const *gt_field(get_variant_format(var).gt);
			for (auto const &[path_idx, sample_idx] : ranges::view::enumerate(psv.representatives_by_path()))
			{
				libbio_assert_neq(SAMPLE_NUMBER_MAX, sample_idx);
				auto const [donor_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(sample_idx);
				auto const &sample(var.samples()[donor_idx]);
				auto const &gt((*gt_field)(sample));
				auto const alt_idx(gt[chr_idx].alt);
				
				// Check whether the ALT was handled. Otherwise don’t modify the zero stored in dst_path_edges.
				if (alt_idx && 0 == unhandled_alt_csum[alt_idx] - unhandled_alt_csum[alt_idx - 1])
				{
					auto const fixed_alt_idx(alt_idx - unhandled_alt_csum[alt_idx]);
					dst_path_edges(var_idx, path_idx) |= fixed_alt_idx;
				}
			}
		}
		
		m_subgraph_variants.clear();
	}
	
	
	void variant_graph_generator::calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx)
	{
		auto const &ref_positions(m_graph.ref_positions());					// 1-based.
		auto &aligned_ref_positions(m_graph.aligned_ref_positions());		// 1-based.
		libbio_assert_lt(1 + node_idx, ref_positions.size());
		libbio_assert_lt(1 + node_idx, aligned_ref_positions.size());
		auto const prev_aligned_pos(aligned_ref_positions[node_idx]);		// Aligned REF position of the /previous/ node.
		auto const prev_ref_pos(ref_positions[node_idx]);					// REF position of the /previous/ node.
		auto const ref_pos(ref_positions[1 + node_idx]);					// REF position of the /given/ node.
		auto const weight(ref_pos - prev_ref_pos);
		aligned_ref_positions[1 + node_idx] = prev_aligned_pos + weight;	// Aligned REF position of the /given/ node.
	}
	
	
	// Take the in-side-edge aligned positions in account.
	void variant_graph_generator::calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx, std::size_t const max_in_alt_edge_aligned_pos)
	{
		calculate_aligned_ref_pos_for_new_node(node_idx);
		
		auto &aligned_ref_positions(m_graph.aligned_ref_positions());	 // 1-based.
		aligned_ref_positions[1 + node_idx] = std::max(aligned_ref_positions[1 + node_idx], max_in_alt_edge_aligned_pos);
	}
	
	
	void variant_graph_generator::assign_alt_edge_labels_and_queue(libbio::variant const &var, std::size_t const node_idx, std::size_t const alt_edge_start_idx)
	{
		auto &alt_edge_labels(m_graph.alt_edge_labels());
		auto const &edge_count_csum(m_graph.alt_edge_count_csum());
		libbio_assert_lt(1 + node_idx, edge_count_csum.size());
		auto const alt_edge_end_idx(edge_count_csum[1 + node_idx]);						// 1-based.
		auto const start_aligned_pos(m_graph.aligned_ref_positions()[1 + node_idx]);	// 1-based.
		
		std::size_t max_alt_edge_aligned_dst_pos(0);
		std::size_t i(alt_edge_start_idx);
		for (auto const &alt : var.alts())
		{
			if (can_handle_variant_alt(alt))
			{
				libbio_assert_lt(i, alt_edge_end_idx);
				alt_edge_labels[i] = alt.alt;
				max_alt_edge_aligned_dst_pos = std::max(max_alt_edge_aligned_dst_pos, start_aligned_pos + alt.alt.size());
				++i;
			}
		}
		
		m_overlap_stack.emplace(var, node_idx, max_alt_edge_aligned_dst_pos);
	}
	
	
	void variant_graph_generator::finalize_graph()
	{
		auto const ref_size(m_reference->size());
		auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(ref_size, 0);
		if (did_create)
			calculate_aligned_ref_pos_for_new_node(node_idx);
	}
}
