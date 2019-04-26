/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <range/v3/all.hpp>
#include <vcf2multialign/can_handle_variant_alts.hh>
#include <vcf2multialign/preprocess/sample_sorter.hh>
#include <vcf2multialign/variant_format.hh>


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	bool sample_sorter::check_variant_for_sample_and_update_state(std::size_t const sample_idx, lb::variant const &var, std::uint8_t const alt_idx)
	{
		// Check the output position of the given sample. REF can always be handled.
		// If the current sample has zero, do not change its end position.
		libbio_assert_neq(0, alt_idx);
		
		auto const &alt(var.alts()[alt_idx - 1]);
		if (can_handle_variant_alt(alt))
		{
			auto const pos(var.zero_based_pos());
			if (m_end_positions_by_sample[sample_idx] <= pos)
			{
				auto const end_pos(variant_end_pos(var));
				libbio_assert_lt(sample_idx, m_end_positions_by_sample.size());
				m_end_positions_by_sample[sample_idx] = end_pos;
				return true;
			}
			else
			{
				// FIXME: logging should be elsewhere.
				std::cerr << "Line " << var.lineno() << ": Sample " << sample_idx << "’s end position was ";
				std::cerr << m_end_positions_by_sample[sample_idx];
				std::cerr << " but the variant’s position was " << pos << ".\n";
			}
		}
		else
		{
			// FIXME: logging should be elsewhere.
			std::cerr << "Line " << var.lineno() << ": Unable to handle ALT " << alt_idx << " (type " << lb::to_string(alt.alt_sv_type) << ").\n";
		}
	
		return false;
	}
	
	
	// Called before processing the next subgraph.
	void sample_sorter::prepare_for_next_subgraph()
	{
		m_path_counter = 1;
		std::fill(m_src_paths.begin(), m_src_paths.end(), 0);
		std::fill(m_end_positions_by_sample.begin(), m_end_positions_by_sample.end(), 0);
	}
	
	
	void sample_sorter::sort_by_variant_and_alt(lb::variant const &var, std::uint8_t const expected_alt_idx)
	{
		libbio_assert_neq(0, expected_alt_idx);

		// FIXME: currently the variable naming is slightly confusing as alt_idx and expected_alt_idx are 1-based but var.alts() expects a zero-based index.
		auto const *gt_field(get_variant_format(var).gt);
		
		// Reset variant-specific state.
		for (auto &word : m_branches_by_path_index.word_range() | ranges::view::take(1 + m_path_counter / branch_vector::ELEMENT_COUNT))
			word.store(0, std::memory_order_release);
	
		// Check which paths branch by setting a bit mask in m_branches_by_path_index.
		// The branching ones will have 0x3 in the end.
		{
			std::size_t sample_idx(0);
			for (auto const path_idx : m_src_paths)
			{
				// Convert to donor and chr indices.
				// If the ALT cannot be handled, use REF in alt_idx_.
				auto const [donor_idx, chr_idx] = m_sample_indexer->donor_and_chr_idx(sample_idx);
				auto const &sample(var.samples()[donor_idx]);
				auto const &gt((*gt_field)(sample));
				auto const alt_idx(gt[chr_idx].alt);
				if (expected_alt_idx == alt_idx)
				{
					auto const alt_idx_(check_variant_for_sample_and_update_state(sample_idx, var, alt_idx) ? alt_idx : 0);
					auto const val((expected_alt_idx == alt_idx_) ? 0x2 : 0x1);
					//std::cerr << "1 donor_idx: " << donor_idx << " chr_idx: " << +chr_idx << " expected_alt_idx: " << +expected_alt_idx << " alt_idx: " << alt_idx << " alt_idx_: " << alt_idx_ << " path_idx: " << path_idx << " val: " << val << std::endl;
					m_branches_by_path_index.fetch_or(path_idx, val, std::memory_order_release);
				}
				else
				{
					m_branches_by_path_index.fetch_or(path_idx, 0x1, std::memory_order_release);
				}

				++sample_idx;
			}
		}
	
		// Assign numbers to the branching paths.
		for (auto const &tup : ranges::view::enumerate(m_branches_by_path_index) | ranges::view::take(m_path_counter))
		{
			auto const [path_idx, val] = tup;
			auto const mask(val.load(std::memory_order_acquire));
			//std::cerr << "2 path_idx: " << path_idx << " mask: " << mask << std::endl;
			if (0x3 == mask)
			{
				libbio_assert_lt(path_idx, m_branching_paths.size());
				m_branching_paths[path_idx] = m_path_counter++;
			}
		}
	
		// Assign path numbers to samples.
		{
			std::size_t sample_idx(0);
			for (auto const src_path_idx : m_src_paths)
			{
				auto const [donor_idx, chr_idx] = m_sample_indexer->donor_and_chr_idx(sample_idx);
				auto const &sample(var.samples()[donor_idx]);
				auto const &gt((*gt_field)(sample));
				auto const alt_idx(gt[chr_idx].alt);
			
				libbio_assert_lt(src_path_idx, m_branching_paths.size());
				auto const dst_path_idx(
					(m_branches_by_path_index.load(src_path_idx, std::memory_order_acquire) == 0x3 && expected_alt_idx == alt_idx)
					? m_branching_paths[src_path_idx]
					: src_path_idx
				);
				libbio_assert_lt(sample_idx, m_dst_paths.size());
				//std::cerr << "3 lineno: " << var.lineno() <<  " sample_idx: " << sample_idx << " dst_path_idx: " << dst_path_idx << std::endl;
				m_dst_paths[sample_idx] = dst_path_idx;
				++sample_idx;
			}
		}
	
		using std::swap;
		swap(m_src_paths, m_dst_paths);
	}
}
