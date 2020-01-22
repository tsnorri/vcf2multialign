/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <range/v3/all.hpp>
#include "algorithms.hh"
#include "overlap_counter.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {

	// Update the overlap count list.
	void overlap_counter::push_count(variant_record const &var, std::size_t const expected_ploidy)
	{
		// By adding the starting and ending position of the new variant,
		// we will end up with two sorted lists the second one having two elements.
		// These can then be merged.

		auto const overlap_count_size(m_overlap_counts.size());
		libbio_assert_lt(0, overlap_count_size);
		
		auto const var_pos(var.variant.zero_based_pos());
		auto const var_end(var_pos + var.size);

		// Calculate the overlap count.
		auto const [var_overlap_count, ploidy] = (count_set_genotype_values(var.variant, 0));
		libbio_always_assert_eq_msg(ploidy, expected_ploidy, "Line ", var.variant.lineno(), ": expected the sample ploidy to match the passed value, got ", ploidy, '.');

		// Don’t add counts if no GT values were set.
		if (0 == var_overlap_count)
			return;
		
		m_overlap_counts.emplace_back(var_pos, var_overlap_count);
		m_overlap_counts.emplace_back(var_end, -1 * var_overlap_count);

		// Sort.
		// Make sure that no counts were added before the front element that is supposed to store the previously calculated running sum.
		libbio_assert_lte(m_overlap_counts.front().position, m_overlap_counts[overlap_count_size].position);
		std::inplace_merge(
			m_overlap_counts.begin(),
			m_overlap_counts.begin() + overlap_count_size,
			m_overlap_counts.end(),
			[](auto const &lhs, auto const &rhs){
				return lhs.position < rhs.position;
			}
		);
		libbio_assert(
			std::is_sorted(
				m_overlap_counts.begin(),
				m_overlap_counts.end(),
				[](auto const &lhs, auto const &rhs){ return lhs.position < rhs.position; }
			)
		);
	}


	void overlap_counter::update_running_sums()
	{
		// Find the ranges of equivalent positions.
		{
			// Don’t join the first overlap_count b.c. the first one was already handled in the previous call to update_overlap_running_sums().
			libbio_assert_neq(m_overlap_counts.begin(), m_overlap_counts.end());
			auto it(m_overlap_counts.begin() + 1);
			while (true)
			{
				auto const res(multiple_adjacent_find(it, m_overlap_counts.end(), [](auto const &oc){
					return oc.position;
				}));
				
				// Check for an empty range.
				if (res.first == res.second)
					break;
				it = res.first;
				auto const end(res.second);
				
				// Removing the items in the end of this loop will cause iterators to be invalidated.
				// Hence, calculate the range start position for later reference.
				auto const range_start_idx(std::distance(m_overlap_counts.begin(), it));
				auto const overlap_count(std::accumulate(it, end, std::int32_t(0), [](auto const sum, auto const &oc) -> std::int32_t {
					return sum + oc.count;
				}));
				
				// Update the count.
				it->count = overlap_count;
				
				// Remove unneeded items and update the range start iterator.
				m_overlap_counts.erase(++it, end);
				it = m_overlap_counts.begin() + range_start_idx + 1;
			}
		}
		
		// Update the sums.
		libbio_assert(!m_overlap_counts.empty());
		std::int32_t current_sum(m_overlap_counts.front().running_sum);
		for (auto &oc : m_overlap_counts | rsv::tail)
		{
			current_sum += oc.count;
			oc.running_sum = current_sum;
			libbio_assert_lte(0, current_sum);
		}
	}
	
	
	void overlap_counter::clean_up_counts(std::size_t const first_unhandled_pos)
	{
		libbio_assert(!m_overlap_counts.empty());
		auto const begin(m_overlap_counts.begin());
		auto const end(m_overlap_counts.end());
		auto const it(std::partition_point(begin, end, [first_unhandled_pos](auto const &oc){
			return oc.position < first_unhandled_pos;
		}));
		if (begin != it)
			m_overlap_counts.erase(begin, it - 1);
		libbio_assert(!m_overlap_counts.empty());
	}
	
	
	auto overlap_counter::max_overlaps_in_range(
		const_iterator overlap_it,
		const_iterator const overlap_end,
		std::size_t const var_pos,
		std::size_t const var_end_pos
	) const -> std::pair <std::int32_t, const_iterator> 
	{
		struct {
			typedef std::pair <std::size_t, std::size_t> interval;
			// We would like to skip the oc that is exactly at the segment start.
			bool operator()(overlap_count const &oc, interval const &ival) const { return oc.position < ival.first; }
			bool operator()(interval const &ival, overlap_count const &oc) const { return ival.second <= oc.position; }
		} overlap_cmp;
		auto const overlap_it_pair(
			std::equal_range(
				overlap_it,
				overlap_end,
				std::make_pair(var_pos, var_end_pos),
				overlap_cmp
			)
		);
		auto const max_overlaps_in_range(
			(overlap_end == overlap_it_pair.first || overlap_it_pair.first == overlap_it_pair.second)
			? std::int32_t(0)
			: ranges::max(ranges::subrange(overlap_it_pair.first, overlap_it_pair.second) | rsv::transform([](auto const &oc){ return oc.running_sum; }))
		);
		auto const max_overlaps(std::max(overlap_it == m_overlap_counts.begin() ? 0 : (overlap_it - 1)->running_sum, max_overlaps_in_range));
		return std::pair <std::int32_t, const_iterator>(max_overlaps, overlap_it_pair.first);
	}
}
