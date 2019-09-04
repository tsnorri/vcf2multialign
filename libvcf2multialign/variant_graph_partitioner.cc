/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/dispatch.hh>
#include <libbio/utility.hh>
#include <vcf2multialign/preprocess/variant_graph_partitioner.hh>
#include <vcf2multialign/utility/can_handle_variant_alts.hh>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	bool variant_graph_partitioner::partition(
		std::vector <std::string> const &field_names_for_filter_by_assigned,
		cut_position_list &out_cut_positions
	)
	{
		m_reader->reset();
		m_reader->set_parsed_fields(lb::vcf_field::ALL);
		
		// Get the field descriptors needed for accessing the values.
		auto const *end_field(m_reader->get_end_field_ptr());
		
		// Determine the fields used for filtering.
		std::vector <lb::vcf_info_field_base *> filter_by_assigned;
		{
			auto const &fields(m_reader->info_fields());
			for (auto const &name : field_names_for_filter_by_assigned)
			{
				auto const it(fields.find(name));
				if (fields.end() == it)
				{
					m_delegate->variant_processor_no_field_for_identifier(name);
					continue;
				}
				
				filter_by_assigned.emplace_back(it->second.get());
			}
		}
		
		std::vector <cut_position> cut_position_tree;
		std::list <dp_ctx> closable_partitions, unclosable_partitions;
		
		{
			// Set up the first cut position and segment.
			cut_position_tree.emplace_back();
			auto &ctx(unclosable_partitions.emplace_back(*m_reader, m_sample_indexer));
			ctx.start_position_idx = cut_position_tree.size() - 1;
		}
		
		auto &handled_line_numbers(out_cut_positions.handled_line_numbers);
		handled_line_numbers.clear();
		
		bool should_continue(false);
		std::size_t overlap_end{};
		do {
			m_reader->fill_buffer();
			should_continue = m_reader->parse(
				[
					this,
					end_field,
					&filter_by_assigned,
					&cut_position_tree,
					&closable_partitions,
					&unclosable_partitions,
					&handled_line_numbers,
				 	&overlap_end
				](lb::transient_variant const &var) -> bool
				{
					auto const lineno(var.lineno());
					auto const var_pos(var.zero_based_pos());
					
					// Check the chromosome name.
					if (var.chrom_id() != m_chromosome_name)
						goto end;
					
					if (! (var_pos < m_reference->size()))
					{
						m_delegate->variant_processor_found_variant_with_position_greater_than_reference_length(var);
						return false;
					}
					
					if (!can_handle_variant_alts(var))
					{
						m_delegate->variant_processor_found_variant_with_no_suitable_alts(var);
						goto end;
					}
					
					// Filter.
					for (auto const *field_ptr : filter_by_assigned)
					{
						if (field_ptr->has_value(var))
						{
							m_delegate->variant_processor_found_filtered_variant(var, *field_ptr);
							goto end;
						}
					}
					
					// Compare the REF column against the reference sequence.
					{
						auto const &ref_col(var.ref());
						std::string_view const ref_sub(m_reference->data() + var_pos, ref_col.size());
						if (ref_col != ref_sub)
						{
							m_delegate->variant_processor_found_variant_with_ref_mismatch(var, ref_sub);
							goto end;
						}
					}
					
					// Variant passes the checks, handle it.
					handled_line_numbers.emplace_back(lineno);
					
					// Check if the current node is a candidate for splitting.
					if (overlap_end <= var_pos)
					{
						check_closable(var_pos, cut_position_tree, unclosable_partitions, closable_partitions);
						
						// If there is a closable partition, make a copy of it.
						if (!closable_partitions.empty())
						{
							auto &new_ctx(unclosable_partitions.emplace_back(*m_reader, m_sample_indexer));
							new_ctx.chain_previous(closable_partitions.front(), var_pos, cut_position_tree);
						}
					}
					
					// Update the scores.
					{
						auto const alt_count(var.alts().size());
						for (auto &ctx : closable_partitions)
							ctx.count_paths(var, alt_count);
						
						lb::parallel_for_each_range_view(
							unclosable_partitions,
							8,
							[&var, alt_count](auto &ctx, std::size_t const j)
							{
								ctx.count_paths(var, alt_count);
							}
						);
					}
					
					{
						// FIXME: save the processed variants’ identifiers / line numbers / something.
						auto const var_end(lb::variant_end_pos(var, *end_field));
						overlap_end = std::max(overlap_end, var_end);
					}
					
				end:
					++m_processed_count;
					if (0 == m_processed_count % 100)
						std::cerr << m_processed_count << '\n';
					return true;
				}
			);
		} while (should_continue);
		
		auto const ref_len(m_reference->size());
		check_closable(ref_len, cut_position_tree, unclosable_partitions, closable_partitions);
		
		// closable_partitions should now contain at most one partitioning.
		if (closable_partitions.empty())
			return false;
		
		{
			auto &positions(out_cut_positions.positions);
			positions.clear();
			auto &ctx(closable_partitions.front());
			auto idx(ctx.start_position_idx);
			while (SIZE_MAX != idx)
			{
				auto const &cut_pos(cut_position_tree[idx]);
				positions.emplace_back(cut_pos.value);
				idx = cut_pos.previous_idx;
			}
			std::reverse(positions.begin(), positions.end());
			out_cut_positions.max_segment_size = ctx.max_size;
			return true;
		}
	}
	
	
	void variant_graph_partitioner::check_closable(
		std::size_t const var_pos,
		std::vector <cut_position> const &cut_position_tree,
		std::list <dp_ctx> &unclosable_partitions,
		std::list <dp_ctx> &closable_partitions
	)
	{
		// Check if there are partitions that can be closed.
		{
			auto it(unclosable_partitions.begin());
			auto const end(unclosable_partitions.end());
			while (it != end)
			{
				auto const &ctx(*it);
				auto const idx(ctx.start_position_idx);
				auto const &cut_pos(cut_position_tree[idx]);
				if (cut_pos.value + m_minimum_subgraph_distance <= var_pos)
				{
					auto const src(it++);
					closable_partitions.splice(closable_partitions.end(), unclosable_partitions, src);
				}
				else
				{
					++it;
				}
			}
		}
		
		// Check if there are more than one closable partitions.
		// There should be no more than one.
		if (1 < closable_partitions.size())
		{
			// Find the smallest partitions.
			path_number_type smallest_size(PATH_NUMBER_MAX);
			for (auto const &ctx : closable_partitions)
				smallest_size = std::min(smallest_size, ctx.max_size);
			
			std::list <dp_ctx> helper;
			auto it(closable_partitions.begin());
			auto const end(closable_partitions.end());
			while (it != end)
			{
				auto &ctx(*it);
				if (ctx.max_size == smallest_size)
				{
					auto const src(it++);
					helper.splice(helper.end(), closable_partitions, src);
				}
				else
				{
					++it;
				}
			}
			
			// Sort by cut position count, keep the most recently cut. This does not worsen the results since
			// we already know that the partition is closable.
			helper.sort([&cut_position_tree](auto const &lhs, auto const &rhs){
				return cut_position_tree[lhs.start_position_idx].value > cut_position_tree[rhs.start_position_idx].value;
			});
			closable_partitions.clear();
			closable_partitions.splice(closable_partitions.end(), helper, helper.begin());
		}
	}
	
	
	std::ostream &operator<<(std::ostream &os, variant_graph_partitioner::dp_ctx const &ctx)
	{
		os << "dp_ctx max_size: " << ctx.max_size << " start_position_idx: " << ctx.start_position_idx;
		return os;
	}
	
	
	void variant_graph_partitioner::dp_ctx::output_reversed_path(std::ostream &os, std::vector <cut_position> const &cut_position_tree) const
	{
		os << "dp_ctx max_size: " << max_size << " path:";
		auto idx(start_position_idx);
		while (SIZE_MAX != idx)
		{
			auto const &cut_pos(cut_position_tree[idx]);
			os << ' ' << cut_pos.value;
			idx = cut_pos.previous_idx;
		}
	}
}
