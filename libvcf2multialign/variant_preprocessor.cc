/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/algorithm.hh>
#include <libbio/bits.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/transform.hpp>
#include <vcf2multialign/preprocess/variant_preprocessor.hh>
#include <vcf2multialign/utility/can_handle_variant_alts.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void variant_preprocessor::update_sample_names()
	{
		// XXX I don’t remember why vcf_reader outputs the sample names as a map with names as keys and positions as indices but it should be fixed unless a good reason is found not to.
		auto const &sample_name_map(m_reader->sample_names());
		m_sample_names.resize(sample_name_map.size());
		for (auto const & [sample_name, idx1] : sample_name_map)
			m_sample_names[idx1 - 1] = sample_name; // Copy.
	}
	
	
	void variant_preprocessor::process(std::vector <std::string> const &field_names_for_filter_by_assigned)
	{
		update_sample_names();
		m_graph.clear();
		
		m_processed_count = 0;
		m_reader->reset();
		m_reader->set_parsed_fields(lb::vcf_field::ALL);
		std::size_t overlap_start(0);
		std::size_t overlap_end(0);
		std::size_t prev_overlap_end(0);
		
		// Get the field descriptors needed for accessing the values.
		libbio::vcf_info_field_end const *end_field{};
		m_reader->get_info_field_ptr("END",	end_field);
		
		// Determine the fields used for filtering.
		std::vector <lb::vcf_info_field_base *> filter_by_assigned;
		{
			auto const &fields(m_reader->info_fields());
			for (auto const &name : field_names_for_filter_by_assigned)
			{
				auto const it(fields.find(name));
				if (fields.end() == it)
				{
					m_delegate->variant_preprocessor_no_field_for_identifier(name);
					continue;
				}
				
				filter_by_assigned.emplace_back(it->second.get());
			}
		}
		
		// Process the variants.
		bool should_continue(false);
		do {
			m_reader->fill_buffer();
			should_continue = m_reader->parse(
				[
					this,
					end_field,
					&filter_by_assigned,
				 	&overlap_start,
				 	&overlap_end,
				 	&prev_overlap_end
				](lb::transient_variant const &var) -> bool
				{
					auto const lineno(var.lineno());
					auto const var_pos(var.zero_based_pos());

					// Check the chromosome name.
					if (var.chrom_id() != m_chromosome_name)
						goto end;
					
					if (! (var_pos < m_reference->size()))
					{
						m_delegate->variant_preprocessor_found_variant_with_position_greater_than_reference_length(var);
						return false;
					}

					if (!can_handle_variant_alts(var))
					{
						m_delegate->variant_preprocessor_found_variant_with_no_suitable_alts(var);
						goto end;
					}
					
					// Filter.
					for (auto const *field_ptr : filter_by_assigned)
					{
						if (field_ptr->has_value(var))
						{
							m_delegate->variant_preprocessor_found_filtered_variant(var, *field_ptr);
							goto end;
						}
					}
					
					// Compare the REF column against the reference sequence.
					{
						auto const &ref_col(var.ref());
						std::string_view const ref_sub(m_reference->data() + var_pos, ref_col.size());
						if (ref_col != ref_sub)
						{
							m_delegate->variant_preprocessor_found_variant_with_ref_mismatch(var, ref_sub);
							goto end;
						}
					}
					
					// Variant passes the checks, handle it.
					{
						auto const var_end(lb::variant_end_pos(var, *end_field));
						
						if (overlap_end + m_minimum_subgraph_distance <= var_pos)
						{
							process_subgraph(prev_overlap_end);
							overlap_start = var_pos;
							prev_overlap_end = overlap_end;
						}
						
						overlap_end = std::max(overlap_end, var_end);
						m_subgraph_variants.emplace_back(var);
					}
					
				end:
					++m_processed_count;
					return true;
				}
			);
		} while (should_continue);
		
		// Process the remaining variants.
		process_subgraph(prev_overlap_end);
		
		// Copy the sample names and finalize.
		m_graph.sample_names() = m_sample_names;
		finalize_graph();
	}
	
	
	// Retrieve the accumulated group of variants and pass them to the worker thread for processing.
	void variant_preprocessor::process_subgraph(std::size_t const prev_overlap_end)
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
		m_delegate->variant_preprocessor_will_handle_subgraph(m_subgraph_variants.size(), path_count);
		
		// Begin a new subgraph; place the beginning between the last variant of the previous subgraph and the first variant of the next subgraph.
		std::size_t const subgraph_start_pos(prev_overlap_end + std::ceil((m_subgraph_variants.front().zero_based_pos() - prev_overlap_end) / 2.0));
		auto const [subgraph_start_node_idx, subgraph_start_alt_edge_start_idx, did_create_node_for_subgraph_start] = m_graph.add_main_node(subgraph_start_pos, 0);
		auto const subgraph_idx(m_graph.add_subgraph(subgraph_start_node_idx, sample_count, m_subgraph_variants.size(), path_count));
		if (did_create_node_for_subgraph_start)
			calculate_aligned_ref_pos_for_new_node(subgraph_start_node_idx);
		
		// Iterate over the variants, collect almost everything except for the ALT indices for the paths.
		libbio_always_assert(m_overlap_stack.empty());
		for (auto const &var : m_subgraph_variants)
		{
			auto const ref_pos(var.zero_based_pos());
			
			// Calculate new indices for the handled ALTs.
			auto const total_alt_count(var.alts().size());
			m_unhandled_alt_csum.clear();
			m_unhandled_alt_csum.resize(1 + total_alt_count, 0);
			for (auto const &[idx, alt] : ranges::view::enumerate(var.alts()))
			{
				auto const prev(m_unhandled_alt_csum[idx]);
				if (can_handle_variant_alt(alt))
					m_unhandled_alt_csum[1 + idx] = prev;
				else
					m_unhandled_alt_csum[1 + idx] = 1 + prev;
			}
			
			auto const handled_alt_count(total_alt_count - m_unhandled_alt_csum.back());
			
			// If there are no ALT edges that would possibly end before the new node (resulting from handing the
			// current variant), create a new node.
			if (m_overlap_stack.empty())
			{
				auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(ref_pos, handled_alt_count);
				if (did_create)
					calculate_aligned_ref_pos_for_new_node(node_idx);
				assign_alt_edge_labels_and_queue(var, node_idx, alt_edge_start_idx);
			}
			else
			{
				// Otherwise check first if the ALT edges need to be handled.
				do
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
						break;
					}
					else if (ref_pos == prev_end_pos)
					{
						auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(prev_end_pos, handled_alt_count);
						if (did_create)
							calculate_aligned_ref_pos_for_new_node(node_idx, top_entry.max_alt_edge_aligned_dst_pos);
						m_graph.connect_alt_edges(top_entry.node_number, node_idx);
						assign_alt_edge_labels_and_queue(var, node_idx, alt_edge_start_idx);
						m_overlap_stack.pop();
						break;
					}
					else // ref_pos > prev_end_pos
					{
						auto const [prev_node_idx, prev_alt_edge_start_idx, did_create] = m_graph.add_main_node(prev_end_pos, 0);
						if (did_create)
							calculate_aligned_ref_pos_for_new_node(prev_node_idx, top_entry.max_alt_edge_aligned_dst_pos);
						m_graph.connect_alt_edges(top_entry.node_number, prev_node_idx);
						m_overlap_stack.pop();
					}
				} while (!m_overlap_stack.empty());
			}
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
		for (auto const &[path_idx, sample_idx] : ranges::view::enumerate(psv.representatives_by_path()))
		{
			libbio_assert_neq(SAMPLE_NUMBER_MAX, sample_idx);
			auto const [donor_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(sample_idx);
			for (auto const &[var_idx, var] : ranges::view::enumerate(m_subgraph_variants))
			{
				auto const *gt_field(get_variant_format(var).gt);
				auto const &sample(var.samples()[donor_idx]);
				auto const &gt((*gt_field)(sample));
				auto const alt_idx(gt[chr_idx].alt);
				
				// Check whether the ALT was handled. Otherwise don’t modify the zero stored in dst_path_edges.
				if (alt_idx && m_unhandled_alt_csum[alt_idx - 1] == m_unhandled_alt_csum[alt_idx])
				{
					auto const fixed_alt_idx(alt_idx - m_unhandled_alt_csum[alt_idx]); // 0 for alt_idx stands for REF which is set to zero in the csum anyway.
					dst_path_edges(path_idx, var_idx) |= fixed_alt_idx;
				}
			}
		}
		
		m_subgraph_variants.clear();
	}
	
	
	void variant_preprocessor::calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx)
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
	void variant_preprocessor::calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx, std::size_t const max_in_alt_edge_aligned_pos)
	{
		calculate_aligned_ref_pos_for_new_node(node_idx);
		
		auto &aligned_ref_positions(m_graph.aligned_ref_positions());	 // 1-based.
		aligned_ref_positions[1 + node_idx] = std::max(aligned_ref_positions[1 + node_idx], max_in_alt_edge_aligned_pos);
	}
	
	
	void variant_preprocessor::assign_alt_edge_labels_and_queue(libbio::variant const &var, std::size_t const node_idx, std::size_t const alt_edge_start_idx)
	{
		auto &alt_edge_labels(m_graph.alt_edge_labels());
		auto const &edge_count_csum(m_graph.alt_edge_count_csum());
		libbio_assert_lt(1 + node_idx, edge_count_csum.size());
		auto const alt_edge_end_idx(edge_count_csum[1 + node_idx]);					// 1-based.
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
	
	
	void variant_preprocessor::finalize_graph()
	{
		auto const ref_size(m_reference->size());
		auto const [node_idx, alt_edge_start_idx, did_create] = m_graph.add_main_node(ref_size, 0);
		if (did_create)
			calculate_aligned_ref_pos_for_new_node(node_idx);
	}
}
