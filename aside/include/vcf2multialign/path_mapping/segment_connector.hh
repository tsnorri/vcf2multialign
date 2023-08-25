/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PATH_MAPPING_SEGMENT_CONNECTOR_HH
#define VCF2MULTIALIGN_PATH_MAPPING_SEGMENT_CONNECTOR_HH

#include <set>
#include <vcf2multialign/path_mapping/types.hh>


namespace vcf2multialign { namespace path_mapping {
	
	class segment_connector
	{
	protected:
		substring_index_multiset	m_substrings_available_lhs;	// Substrings available in the lhs segment.
		std::size_t					m_founder_count{};
		std::size_t					m_slots_available_lhs{};
		
	public:
		segment_connector(std::size_t const founder_count):
			m_founder_count(founder_count)
		{
		}
		
		void setup_with_complete_segment(std::size_t const initial_lhs_count);
		void setup_with_selected_substrings(substring_index_multiset &&substrings_available_lhs);
		
		void make_edges(
			path_item_vector const &path_counts,
			std::size_t const substring_count_rhs,
			edge_vector &edges,
			substring_index_vector &substrings_added_to_lhs
		);
		
		void make_edges_using_specific_substrings(
			path_item_vector const &path_counts,
			substring_index_multiset &&substrings_available,
			edge_vector &edges,
			substring_index_vector &substrings_added_to_lhs
		);
	};
}}

#endif