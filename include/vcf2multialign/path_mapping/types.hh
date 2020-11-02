/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PATH_MAPPING_TYPES_HH
#define VCF2MULTIALIGN_PATH_MAPPING_TYPES_HH

#include <ostream>
#include <vector>


namespace vcf2multialign { namespace path_mapping {
	
	typedef std::uint32_t							substring_index_type;
	typedef substring_index_type					substring_count_type;
	typedef std::vector <substring_index_type>		substring_index_vector;
	
	enum { UNASSIGNED_INDEX = std::numeric_limits <substring_index_type>::max() };
	
	struct edge; // Fwd.
	typedef std::vector <edge>		edge_vector;
	
	struct path_item; // Fwd.
	typedef std::vector <path_item>	path_item_vector;
	
	
	struct edge
	{
		substring_index_type lhs_idx{};
		substring_index_type rhs_idx{};
		
		edge() = default;
		
		edge(substring_index_type lhs_idx_, substring_index_type rhs_idx_):
			lhs_idx(lhs_idx_),
			rhs_idx(rhs_idx_)
		{
		}
		
		auto to_tuple() const { return std::make_tuple(lhs_idx, rhs_idx); }
	};
	
	inline bool operator<(edge const &lhs, edge const &rhs)
	{
		return lhs.to_tuple() < rhs.to_tuple();
	}

	inline std::ostream &operator<<(std::ostream &os, edge const &edge)
	{
		os << '(' << edge.lhs_idx << ", " << edge.rhs_idx << ')';
		return os;
	}
}}

#endif
