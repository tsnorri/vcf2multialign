/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/dispatch.hh>
#include <libbio/utility.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/utility/can_handle_variant_alts.hh>

namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	bool variant_partitioner::partition(
		std::vector <std::string> const &field_names_for_filter_by_assigned,
		preprocessing_result &out_preprocessing_result,
		bool const should_start_from_current_variant
	)
	{
		this->m_reader->set_parsed_fields(vcf::field::ALL);
		
		// Get the field descriptors needed for accessing the values.
		auto const *end_field(this->m_reader->get_end_field_ptr());
		
		// Determine the fields used for filtering.
		std::vector <vcf::info_field_base *> filter_by_assigned;
		this->fill_filter_by_assigned(field_names_for_filter_by_assigned, filter_by_assigned, *m_delegate);
		
		std::vector <cut_position> cut_position_tree;
		std::list <dp_ctx> closable_partitions, unclosable_partitions;
		
		{
			// Set up the first cut position and segment.
			cut_position_tree.emplace_back();
			auto &ctx(unclosable_partitions.emplace_back(*m_delegate, *this->m_reader, m_sample_indexer));
			ctx.start_position_idx = cut_position_tree.size() - 1;
		}
		
		auto &handled_line_numbers(out_preprocessing_result.handled_line_numbers);
		handled_line_numbers.clear();
		
		std::size_t overlap_end{};
		auto handling_callback(
			[
				this,
				end_field,
				&filter_by_assigned,
				&cut_position_tree,
				&closable_partitions,
				&unclosable_partitions,
				&handled_line_numbers,
			 	&overlap_end
			](vcf::transient_variant const &var) -> bool
			{
				auto const lineno(var.lineno());
				auto const var_pos(var.zero_based_pos());
				auto const var_end(vcf::variant_end_pos(var, *end_field));
				libbio_always_assert_lte(var_pos, var_end);
				
				switch (this->check_variant(var, filter_by_assigned, *m_delegate))
				{
					case variant_check_status::PASS:
						break;
					case variant_check_status::ERROR:
						goto end;
					case variant_check_status::FATAL_ERROR:
						return false;
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
						auto &new_ctx(unclosable_partitions.emplace_back(*m_delegate, *this->m_reader, m_sample_indexer));
						new_ctx.chain_previous(closable_partitions.front(), var_pos, cut_position_tree);
					}
				}
				
				// Update the scores.
				{
					// FIXME: count_paths expects a vcf::variant. Otherwise the copy would not be needed.
					vcf::variant persistent_var(var);
					
					auto const alt_count(var.alts().size());
					for (auto &ctx : closable_partitions)
						ctx.count_paths(persistent_var, alt_count);
					
					lb::parallel_for_each_range_view(
						unclosable_partitions,
						8,
						[&persistent_var, alt_count](auto &ctx, std::size_t const j)
						{
							ctx.count_paths(persistent_var, alt_count);
						}
					);
				}
				
				{
					auto const var_end(vcf::variant_end_pos(var, *end_field));
					overlap_end = std::max(overlap_end, var_end);
				}
				
			end:
				m_processed_count.fetch_add(1, std::memory_order_relaxed);
				return true;
			}
		);
		
		if (should_start_from_current_variant)
			handling_callback(this->m_reader->current_variant());
		this->m_reader->parse(handling_callback);
		
		auto const ref_len(m_reference->size());
		check_closable(ref_len, cut_position_tree, unclosable_partitions, closable_partitions);
		
		// closable_partitions should now contain at most one partitioning.
		if (closable_partitions.empty())
			return false;
		
		{
			auto &positions(out_preprocessing_result.positions);
			positions.clear();
			positions.emplace_back(ref_len);
			
			auto &ctx(closable_partitions.front());
			auto idx(ctx.start_position_idx);
			while (SIZE_MAX != idx)
			{
				auto const &cut_pos(cut_position_tree[idx]);
				positions.emplace_back(cut_pos.value);
				idx = cut_pos.previous_idx;
			}
			std::reverse(positions.begin(), positions.end());
			libbio_assert(std::is_sorted(positions.begin(), positions.end()));
			libbio_assert_eq(positions.end(), std::adjacent_find(positions.begin(), positions.end()));
			out_preprocessing_result.max_segment_size = ctx.max_size;
			out_preprocessing_result.is_valid = true; // Mark as valid.
			return true;
		}
	}
	
	
	void variant_partitioner::check_closable(
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
	
	
	std::ostream &operator<<(std::ostream &os, variant_partitioner::dp_ctx const &ctx)
	{
		os << "dp_ctx max_size: " << ctx.max_size << " start_position_idx: " << ctx.start_position_idx;
		return os;
	}
	
	
	void variant_partitioner::dp_ctx::output_reversed_path(std::ostream &os, std::vector <cut_position> const &cut_position_tree) const
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
