/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PATH_MAPPING_PATH_MAPPER_HH
#define VCF2MULTIALIGN_PATH_MAPPING_PATH_MAPPER_HH

#include <set>
#include <vcf2multialign/path_mapping/types.hh>


namespace vcf2multialign { namespace path_mapping {
	
	struct substring_item
	{
		substring_index_type substring_idx{};
		substring_count_type count{};
		
		substring_item() = default;
		
		substring_item(substring_index_type substring_idx_, substring_count_type count_):
			substring_idx(substring_idx_),
			count(count_)
		{
		}
		
		explicit substring_item(substring_index_type substring_idx_):
			substring_item(substring_idx_, 0)
		{
		}
	};
	
	inline bool operator<(substring_item const &lhs, substring_item const &rhs)
	{
		return lhs.substring_idx < rhs.substring_idx;
	}
	
	
	struct path_item
	{
		substring_index_type lhs_idx{};
		substring_index_type rhs_idx{};
		substring_count_type count{};
		
		path_item() = default;
		
		path_item(substring_index_type lhs_idx_, substring_index_type rhs_idx_, substring_count_type count_):
			lhs_idx(lhs_idx_),
			rhs_idx(rhs_idx_),
			count(count_)
		{
		}
		
		path_item(substring_index_type lhs_idx_, substring_index_type rhs_idx_):
			path_item(lhs_idx_, rhs_idx_, 0)
		{
		}
		
		auto to_path_index_tuple() const { return std::make_tuple(lhs_idx, rhs_idx); }
	};
	
	inline bool operator<(path_item const &lhs, path_item const &rhs)
	{
		return lhs.to_path_index_tuple() < rhs.to_path_index_tuple();
	}
	
	
	class path_mapper
	{
	public:
		typedef std::vector <substring_index_vector>	substring_index_inv_mapping;
		typedef std::set <substring_index_type>			substring_index_set;
		
	protected:
		substring_index_inv_mapping	m_idxs_by_substring_lhs{};		// Founder (output stream) indices by substring index.
		substring_index_inv_mapping	m_idxs_by_substring_rhs{};		// Founder (output stream) indices by substring index.
		substring_index_vector		m_idxs_available_lhs{};			// Founder (output stream) indices available.
		substring_index_set			m_idxs_available_rhs{};			// Founder (output stream) indices available.
		substring_index_vector		m_string_idxs_by_founder_lhs{};	// Substring indices by founder index.
#ifndef NDEBUG
		edge_vector					m_prev_edges;
#endif
		std::size_t					m_founder_count{};
		
	public:
		// first_subgraph_size equals the number of paths in the first subgraph.
		path_mapper(std::size_t first_subgraph_size, std::size_t founder_count):
			m_idxs_by_substring_lhs(first_subgraph_size),
			m_idxs_by_substring_rhs(first_subgraph_size),
			m_string_idxs_by_founder_lhs(founder_count),
			m_founder_count(founder_count)
		{
		}
		
		void setup_with_complete_segment(std::size_t const initial_lhs_count);
		void setup_with_selected_substrings(substring_index_multiset const &selected_substrings);
		
		void set_rhs_subgraph_size(std::size_t const count);
		void add_substrings(substring_index_vector const &substrings_added_to_lhs);
		void assign_edges_to_founders(edge_vector const &edges);
		void update_string_indices();
		void end_subgraph();
		
		substring_index_vector const &founder_indices_available_lhs() const { return m_idxs_available_lhs; }
		substring_index_vector const &string_indices_by_founder() const { return m_string_idxs_by_founder_lhs; }
	};
}}

#endif
