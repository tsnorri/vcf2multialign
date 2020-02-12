/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_ALGORITHMS_HH
#define VCF2MULTIALIGN_COMBINE_MSA_ALGORITHMS_HH

#include <range/v3/all.hpp>


namespace vcf2multialign {
	
	// Treat the given range as segments.
	template <typename t_range, typename t_proj_fn, typename t_start_fn, typename t_handle_fn, typename t_end_fn>
	void segment(t_range const &range, t_proj_fn &&proj_fn, t_start_fn &&start_fn, t_handle_fn &&handle_fn, t_end_fn &&end_fn)
	{
		if (range.empty())
			return;
		
		auto current_id(proj_fn(range.front()));
		start_fn(range.front(), current_id);
		handle_fn(range.front(), current_id);
		for (auto const &item : range | ranges::view::tail)
		{
			auto const next_id(proj_fn(item));
			if (current_id != next_id)
			{
				end_fn(current_id);
				start_fn(item, next_id);
			}
			handle_fn(item, next_id);
			current_id = next_id;
		}
		end_fn(current_id);
	}
	
	
	// Find the first subrange of equal consecutive elements in O(n) time.
	// If not found, return an empty range.
	template <typename t_iterator, typename t_proj>
	std::pair <t_iterator, t_iterator> multiple_adjacent_find(t_iterator it, t_iterator end, t_proj &&proj)
	{
		// Check for an empty range.
		if (it == end)
			return {end, end};
		
		auto candidate_segment_begin(it);
		auto candidate_id(proj(*it));
		++it;
		
		// Check for a singular range.
		if (it == end)
			return {end, end};
		
		bool found(false);
		while (it != end)
		{
			auto current_id(proj(*it));
			if (candidate_id == current_id)
			{
				found = true;
				++it;
				continue;
			}
			
			if (found)
				return {candidate_segment_begin, it};
			
			candidate_id = std::move(current_id);
			candidate_segment_begin = it;
			++it;
		}
		
		if (found)
			return {candidate_segment_begin, it};
		return {end, end};
	}
}

#endif
