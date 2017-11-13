/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/graph_range.hh>

namespace vcf2multialign {
	
	std::size_t graph_range::seq_position(std::size_t const var_lineno) const
	{
#ifndef NDEBUG
		auto const indices_size(m_skipped_lines.indices.size());
		auto const rank(m_skipped_lines.index_rank_1_support.rank(indices_size));
		if (indices_size - rank != m_variant_count)
		{
			std::cerr << "indices_size: " << indices_size << " rank: " << rank << " variant_count: " << m_variant_count << std::endl;
			fail();
		}
#endif
		
		// Calculate the index to be used taking skipped variants into account.
		assert(m_start_lineno <= var_lineno);
		std::size_t const line_idx(var_lineno - m_start_lineno);
		assert(line_idx < indices_size);
		always_assert(0 == m_skipped_lines.indices[line_idx]);
		auto const skipped_before_current(m_skipped_lines.index_rank_1_support.rank(line_idx));
		std::size_t const seq_position(line_idx - skipped_before_current);
		return seq_position;
	}
	
	
	bool graph_range::contains_var_lineno(std::size_t const var_lineno) const
	{
		if (var_lineno < m_start_lineno)
			return false;
		
		// The number of elements in m_skipped_lines.indices is equal to
		// the length of this graph range (includeing skipped variants).
		if (var_lineno - m_start_lineno < m_skipped_lines.indices.size())
			return true;
		
		return false;
	}
	
	
	bool graph_range::has_valid_alts(std::size_t const var_lineno) const
	{
		assert(m_start_lineno <= var_lineno);
		std::size_t const line_idx(var_lineno - m_start_lineno);
		assert(line_idx < m_skipped_lines.indices.size());
		return (0 == m_skipped_lines.indices[line_idx]);
	}
}
