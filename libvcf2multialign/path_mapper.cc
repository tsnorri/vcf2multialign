/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/zip.hpp>
#include <vcf2multialign/path_mapping/path_mapper.hh>

namespace lb	= libbio;
namespace rsv	= ranges::view;


namespace vcf2multialign { namespace path_mapping {
	
	void path_mapper::setup(std::size_t const initial_lhs_count)
	{
		libbio_assert_lte(initial_lhs_count, m_idxs_by_substring_lhs.size());
		for (std::size_t i(0); i < initial_lhs_count; ++i)
			m_idxs_by_substring_lhs[i].push_back(i);
		
		// Mark the streams that are not associated with a founder sequence available.
		for (std::size_t i(initial_lhs_count); i < m_founder_count; ++i)
			m_idxs_available_lhs.push_back(i);
	}
	
	
	void path_mapper::add_substrings(substring_index_vector const &substrings_added_to_lhs)
	{
		// Additional substrings (copies) may be assigned to lhs before assign_edges_to_founders is called.
		libbio_assert_lte(substrings_added_to_lhs.size(), m_idxs_available_lhs.size());
		for (auto const substring_idx : substrings_added_to_lhs)
		{
			auto const founder_idx(m_idxs_available_lhs.back());
			m_idxs_available_lhs.pop_back();
			libbio_assert_lt(substring_idx, m_idxs_by_substring_lhs.size());
			m_idxs_by_substring_lhs[substring_idx].push_back(founder_idx);
		}
	}
	
	
	void path_mapper::assign_edges_to_founders(edge_vector const &edges)
	{
		// The edge list has edges of the following types:
		// – Both set to a substring index.
		// – Lhs set to UNASSIGNED_INDEX
		// – Rhs set to UNASSIGNED_INDEX
		// Here we determine the founder indices to which the lhs substrings have been assigned and assign them the rhs substrings.
		// Thus, the precondition is that each founder has an lhs substring. The postcondition is that each founder has an rhs substring, too.
		
		// Update m_idxs_by_substring_lhs s.t. each substring has one or more founder sequences associated with it.
		if (!edges.empty())
		{
			m_idxs_available_rhs.clear();
			for (substring_index_type i(0); i < m_founder_count; ++i)
				m_idxs_available_rhs.insert(i);
			
			auto edge_range_begin(edges.begin());
			while (true)
			{
				// Determine the range of edges with lhs_idx set to current_lhs_substring_idx.
				// Then store the founder indices to m_idxs_by_substring_rhs.
				auto const current_lhs_substring_idx(edge_range_begin->lhs_idx);
				auto const edge_range_end(std::partition_point(edge_range_begin, edges.end(), [current_lhs_substring_idx](auto const &edge){
					return edge.lhs_idx == current_lhs_substring_idx;
				}));
				
				auto const edge_range(ranges::subrange(edge_range_begin, edge_range_end));

				try
				{
					if (UNASSIGNED_INDEX == current_lhs_substring_idx)
					{
						// Lhs index is unassigned. This means that the rhs index can be assigned to any founder.
						// This range should be handled last, so the remaining founders may be used.
						libbio_assert_eq(edges.end(), edge_range_end);
						for (auto const &edge : edge_range)
						{
							libbio_assert_neq(UNASSIGNED_INDEX, edge.rhs_idx);
							libbio_assert(!m_idxs_available_rhs.empty());
							auto const it(m_idxs_available_rhs.begin());
							auto const founder_idx(*it);
							libbio_assert_lt(edge.rhs_idx, m_idxs_by_substring_rhs.size());
							m_idxs_by_substring_rhs[edge.rhs_idx].push_back(founder_idx);
							m_idxs_available_rhs.erase(it);
						}
					}
					else
					{
						libbio_assert_lt(current_lhs_substring_idx, m_idxs_by_substring_lhs.size());
						auto const &founder_idxs_lhs(m_idxs_by_substring_lhs[current_lhs_substring_idx]);
						try
						{
							libbio_assert_lte_msg(
								edge_range.size(),
								founder_idxs_lhs.size(),
								"Expected edge_range.size() to be less than or equal to founder_idxs_lhs.size(), got ",
								edge_range.size(),
								" and ",
								founder_idxs_lhs.size()
							);
						
							for (auto const &[founder_idx, edge] : rsv::zip(founder_idxs_lhs, edge_range))
							{
								if (UNASSIGNED_INDEX != edge.rhs_idx)
								{
									auto const it(m_idxs_available_rhs.find(founder_idx));
									libbio_assert_neq(it, m_idxs_available_rhs.end());
									libbio_assert_lt(edge.rhs_idx, m_idxs_by_substring_rhs.size());
									m_idxs_by_substring_rhs[edge.rhs_idx].push_back(founder_idx);
									m_idxs_available_rhs.erase(it);
								}
							}
						}
						catch (lb::assertion_failure_exception const &)
						{
							std::cerr << "founder_idxs_lhs:";
							for (auto const idx : founder_idxs_lhs)
								std::cerr << ' ' << idx;
							std::cerr << '\n';
							throw;
						}
					}
				}
				catch (lb::assertion_failure_exception const &)
				{
					std::cerr << "current_lhs_substring_idx: " << current_lhs_substring_idx << '\n';
					std::cerr << "edge_range:\n";
					for (auto const &edge : edge_range)
						std::cerr << edge << '\n';
					std::cerr << "edges:\n";
					for (auto const &edge : edges)
						std::cerr << edge << '\n';
#ifndef NDEBUG
					std::cerr << "prev_edges:\n";
					for (auto const &edge : m_prev_edges)
						std::cerr << edge << '\n';
#endif
					std::cerr << "m_idxs_by_substring_lhs:\n";
					for (auto const &[idx, vec] : rsv::enumerate(m_idxs_by_substring_lhs))
					{
						std::cerr << idx << ':';
						for (auto const founder_idx : vec)
							std::cerr << ' ' << founder_idx;
						std::cerr << '\n';
					}

					throw;
				}
				
				if (edges.end() == edge_range_end)
					break;
				
				edge_range_begin = edge_range_end;
			}
		}

#ifndef NDEBUG
		m_prev_edges = edges;
#endif
	}
	
	
	void path_mapper::update_string_indices()
	{
		// m_idxs_available_lhs could be used instead.
		std::fill(m_string_idxs_by_founder_lhs.begin(), m_string_idxs_by_founder_lhs.end(), UNASSIGNED_INDEX);
		
		for (auto const &[substring_idx, founder_idx_vec] : rsv::enumerate(m_idxs_by_substring_lhs))
		{
			for (auto const founder_idx : founder_idx_vec)
			{
				libbio_assert_lt(founder_idx, m_string_idxs_by_founder_lhs.size());
				m_string_idxs_by_founder_lhs[founder_idx] = substring_idx;
			}
		}
	}
	
	
	void path_mapper::end_subgraph()
	{
		using std::swap;
		
		swap(m_idxs_by_substring_lhs, m_idxs_by_substring_rhs);
		
		for (auto &idx_vec : m_idxs_by_substring_rhs)
			idx_vec.clear();
		
		// Copy the available founder indices.
		m_idxs_available_lhs.clear();
		std::copy(m_idxs_available_rhs.begin(), m_idxs_available_rhs.end(), std::back_inserter(m_idxs_available_lhs));
	}
}}
