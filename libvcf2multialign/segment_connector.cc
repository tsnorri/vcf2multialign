/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <range/v3/view/reverse.hpp>
#include <vcf2multialign/path_mapping/path_mapper.hh>
#include <vcf2multialign/path_mapping/segment_connector.hh>

namespace rsv	= ranges::view;


namespace {
	
	template <typename t_collection, typename t_value>
	std::pair <bool, typename t_collection::const_iterator> lower_bound_eq(t_collection const &coll, t_value const &val)
	{
		auto const it(coll.lower_bound(val));
		return {it != coll.end() && *it == val, it};
	}
}


namespace vcf2multialign { namespace path_mapping {
	
	void segment_connector::setup(std::size_t const initial_lhs_count)
	{
		m_slots_available_lhs = m_founder_count - initial_lhs_count;
		
		for (std::size_t i(0); i < initial_lhs_count; ++i)
			m_substrings_available_lhs.insert(i);
	}
	
	
	void segment_connector::make_edges(
		path_item_vector const &path_counts,
		std::size_t const substring_count_rhs,
		edge_vector &edges,
		substring_index_vector &substrings_added_to_lhs
	)
	{
		// Precondition:
		// – m_substrings_available_lhs contains the substring indices that are “available” (to be written to founders) on lhs.
		// – path_counts contains the sorted edge candidates.
		//
		// Postcondition:
		// – m_substrings_available_lhs contains the substring indices that are “available” on lhs on the next iteration, including the indices that were copied.
		// – edges contains the actual edges (no more than there are founders).
		// – substrings_added_to_lhs contains the string indices that were added to lhs on this iteration.
		
		edges.clear();
		substrings_added_to_lhs.clear();
		
		std::multiset <substring_index_type> substrings_available_rhs;	// Substrings available in the rhs segment.
		// Fill with the existing string indices.
		for (substring_index_type i(0); i < substring_count_rhs; ++i)
			substrings_available_rhs.insert(i);
		
		auto new_substrings_available_rhs(substrings_available_rhs); // Copy the existing string indices b.c. they are needed during the next iteration.
		std::size_t slots_available_rhs(m_founder_count - substring_count_rhs); // “Slots” are founders that have not been assigned a substring.
		auto founders_available(m_founder_count);
		
		libbio_assert_lte(m_slots_available_lhs, founders_available);
		libbio_assert_lte(slots_available_rhs, founders_available);
		
		// Create the edges based on path_counts.
		for (auto const &item : rsv::reverse(path_counts))
		{
			auto const [found_lhs, lhs_it] = lower_bound_eq(m_substrings_available_lhs, item.lhs_idx);
			auto const [found_rhs, rhs_it] = lower_bound_eq(substrings_available_rhs,   item.rhs_idx);
			
			if (found_lhs && found_rhs)
			{
				--founders_available;
				m_substrings_available_lhs.erase(lhs_it);	// Erases one.
				substrings_available_rhs.erase(rhs_it);		// Erases one.
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			else if (found_lhs && slots_available_rhs)
			{
				--founders_available;
				m_substrings_available_lhs.erase(lhs_it); // Erases one.
				--slots_available_rhs;
				new_substrings_available_rhs.insert(item.rhs_idx);
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			else if (found_rhs && m_slots_available_lhs)
			{
				--founders_available;
				substrings_available_rhs.erase(rhs_it); // Erases one.
				--m_slots_available_lhs;
				substrings_added_to_lhs.push_back(item.lhs_idx);
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			else if (m_slots_available_lhs && slots_available_rhs)
			{
				--founders_available;
				--m_slots_available_lhs;
				--slots_available_rhs;
				substrings_added_to_lhs.push_back(item.lhs_idx);
				new_substrings_available_rhs.insert(item.rhs_idx);
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			
			// Stop if founder count has been reached.
			if (0 == founders_available)
			{
				libbio_assert_eq(0, m_slots_available_lhs);
				libbio_assert_eq(0, slots_available_rhs);
				goto end;
			}
		}
		
		// Connect the remaining edges arbitrarily while there are available founders.
		while (founders_available && !(m_substrings_available_lhs.empty() || substrings_available_rhs.empty()))
		{
			auto const lhs_begin(m_substrings_available_lhs.begin());
			auto const rhs_begin(substrings_available_rhs.begin());
			edges.emplace_back(*lhs_begin, *rhs_begin);
			m_substrings_available_lhs.erase(lhs_begin);
			substrings_available_rhs.erase(rhs_begin);
			--founders_available;
		}
		
		libbio_assert_lte(m_slots_available_lhs, founders_available);
		libbio_assert_lte(slots_available_rhs, founders_available);
		
		while (m_slots_available_lhs && !m_substrings_available_lhs.empty())
		{
			auto const lhs_begin(m_substrings_available_lhs.begin());
			edges.emplace_back(*lhs_begin, UNASSIGNED_INDEX);
			m_substrings_available_lhs.erase(lhs_begin);
			--m_slots_available_lhs;
		}
		
		while (slots_available_rhs && !substrings_available_rhs.empty())
		{
			auto const rhs_begin(substrings_available_rhs.begin());
			edges.emplace_back(UNASSIGNED_INDEX, *rhs_begin);
			substrings_available_rhs.erase(rhs_begin);
			--slots_available_rhs;
		}
		
	end:
		// Clean up.
		m_slots_available_lhs = slots_available_rhs;
		
		using std::swap;
		swap(m_substrings_available_lhs, new_substrings_available_rhs);

		libbio_assert_lte(edges.size(), m_founder_count);
	}
}}
