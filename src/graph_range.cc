/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/graph_range.hh>

namespace vcf2multialign {
	
	std::size_t graph_range::seq_position(std::size_t const var_lineno) const
	{
		// Calculate the index to be used taking skipped variants into account.
		std::size_t const line_idx(var_lineno - m_start_lineno);
		always_assert(0 == m_skipped_lines.indices[line_idx]);
		auto const skipped_before_current(m_skipped_lines.index_rank_1_support.rank(line_idx));
		std::size_t const seq_position(line_idx - skipped_before_current);
		return seq_position;
	}
	
	
	bool graph_range::contains_var_lineno(std::size_t const var_lineno) const
	{
		if (var_lineno < m_start_lineno)
			return false;
		
		if (var_lineno - m_start_lineno < m_skipped_lines.indices.size())
			return true;
		
		return false;
	}
}
