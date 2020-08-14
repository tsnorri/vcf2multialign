/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */


#include <range/v3/view/subrange.hpp>
#include "founder_sequence_greedy_generator.hh"
#include "utility.hh"

// Substring refers here to a labelled graph path through one subgraph whereas path refers to a graph path through two consecutive subgraphs.

namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef std::vector <bool>						bit_vector;
	
	typedef std::uint32_t							substring_index_type;
	typedef substring_index_type					substring_count_type;
	typedef std::vector <substring_index_type>		substring_index_vector;
	typedef std::set <substring_index_type>			substring_index_set;
	typedef std::vector <substring_index_vector>	substring_index_inv_mapping;
	
	enum { UNASSIGNED_INDEX = std::numeric_limits <substring_index_type>::max() };
	
	struct edge; // Fwd.
	typedef std::vector <edge>		edge_vector;
	
	struct path_item; // Fwd.
	typedef std::vector <path_item>	path_item_vector;
	
	template <typename t_collection, typename t_value>
	std::pair <bool, typename t_collection::const_iterator> lower_bound_eq(t_collection const &coll, t_value const &val)
	{
		auto const it(coll.lower_bound(val));
		return {it != coll.end() && *it == val, it};
	}
	
	
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
	
	bool operator<(edge const &lhs, edge const &rhs)
	{
		return lhs.to_tuple() < rhs.to_tuple();
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
	
	bool operator<(path_item const &lhs, path_item const &rhs)
	{
		return lhs.to_path_index_tuple() < rhs.to_path_index_tuple();
	}
	
	
	class segment_connector
	{
	protected:
		std::multiset <substring_index_type>	m_substrings_available_lhs;	// Substrings available in the lhs segment.
		std::size_t								m_founder_count{};
		std::size_t								m_slots_available_lhs{};
		
	public:
		segment_connector(std::size_t const founder_count):
			m_founder_count(founder_count)
		{
		}
		
		void setup(std::size_t const initial_lhs_count);
		
		void make_edges(
			path_item_vector const &path_counts,
			std::size_t const substring_count_rhs,
			edge_vector &edges,
			substring_index_vector &substrings_added_to_lhs
		);
	};
	
	
	class path_mapper
	{
	protected:
		substring_index_inv_mapping	m_idxs_by_substring_lhs{};		// Founder (output stream) indices by substring index.
		substring_index_inv_mapping	m_idxs_by_substring_rhs{};		// Founder (output stream) indices by substring index.
		substring_index_vector		m_idxs_available_lhs{};			// Founder (output stream) indices available.
		substring_index_set			m_idxs_available_rhs{};			// Founder (output stream) indices available.
		substring_index_vector		m_string_idxs_by_founder_lhs{};	// Substring indices by founder index.
		std::size_t					m_founder_count{};
		
	public:
		path_mapper(std::size_t founder_count):
			m_idxs_by_substring_lhs(founder_count), // Distinct substring count ≤ founder_count.
			m_idxs_by_substring_rhs(founder_count),
			m_string_idxs_by_founder_lhs(founder_count),
			m_founder_count(founder_count)
		{
		}
		
		void setup(std::size_t const initial_lhs_count);
		
		void add_substrings(substring_index_vector const &substrings_added_to_lhs);
		void assign_edges_to_founders(edge_vector const &edges);
		void update_string_indices();
		void end_subgraph();
		
		substring_index_vector const &founder_indices_available_lhs() const { return m_idxs_available_lhs; }
		substring_index_vector const &string_indices_by_founder() const { return m_string_idxs_by_founder_lhs; }
	};
	
	
	class sequence_writer
	{
	protected:
		std::string_view			m_reference;
		v2m::variant_graph const	*m_graph{};
		
	public:
		sequence_writer(v2m::variant_graph const &graph, v2m::vector_type const &reference):
			m_reference(reference.data(), reference.size()),
			m_graph(&graph)
		{
		}
		
		void output_ref(std::size_t const node_idx, std::ostream &os) const { output_ref(node_idx, 1 + node_idx, os); }
		void output_subgraph_path(std::size_t const subgraph_idx, substring_index_type const path_idx, std::ostream &os) const;
		void output_ref_for_subgraph(std::size_t const subgraph_idx, std::ostream &os) const;
		void output_gaps_for_subgraph(std::size_t const subgraph_idx, std::ostream &os) const { output_char_for_subgraph(subgraph_idx, '-', os); }
		void output_n_for_subgraph(std::size_t const subgraph_idx, std::ostream &os) const { output_char_for_subgraph(subgraph_idx, 'N', os); }
		
	protected:
		void output_char_for_subgraph(std::size_t const subgraph_idx, char const c, std::ostream &os) const;
		void output_ref(substring_index_type const node_idx, substring_index_type const next_node_idx, std::ostream &os) const;
		void output_alt(substring_index_type const node_idx, substring_index_type const next_node_idx, std::size_t const alt_idx, std::ostream &os) const;
	};
	
	
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
				m_substrings_available_lhs.erase(lhs_it);	// Erases one.
				substrings_available_rhs.erase(rhs_it);		// Erases one.
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			else if (found_lhs && slots_available_rhs)
			{
				m_substrings_available_lhs.erase(lhs_it); // Erases one.
				--slots_available_rhs;
				new_substrings_available_rhs.insert(item.rhs_idx);
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			else if (found_rhs && m_slots_available_lhs)
			{
				substrings_available_rhs.erase(rhs_it); // Erases one.
				--m_slots_available_lhs;
				substrings_added_to_lhs.push_back(item.lhs_idx);
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			else if (m_slots_available_lhs && slots_available_rhs)
			{
				--m_slots_available_lhs;
				--slots_available_rhs;
				substrings_added_to_lhs.push_back(item.lhs_idx);
				new_substrings_available_rhs.insert(item.rhs_idx);
				edges.emplace_back(item.lhs_idx, item.rhs_idx);
			}
			
			// Stop if founder count has been reached.
			--founders_available;
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
	}
	
	
	void path_mapper::setup(std::size_t const initial_lhs_count)
	{
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
						m_idxs_by_substring_rhs[edge.rhs_idx].push_back(founder_idx);
						m_idxs_available_rhs.erase(it);
					}
				}
				else
				{
					auto const &founder_idxs_lhs(m_idxs_by_substring_lhs[current_lhs_substring_idx]);
					libbio_assert_eq_msg(
						edge_range.size(),
						founder_idxs_lhs.size(),
						"Expected edge_range.size() to be equal to founder_idxs_lhs.size(), got ",
						edge_range.size(),
						" and ",
						founder_idxs_lhs.size(),
						"."
					);
				
					for (auto const &[founder_idx, edge] : rsv::zip(founder_idxs_lhs, edge_range))
					{
						if (UNASSIGNED_INDEX != edge.rhs_idx)
						{
							auto const it(m_idxs_available_rhs.find(founder_idx));
							libbio_assert_neq(it, m_idxs_available_rhs.end());
							m_idxs_by_substring_rhs[edge.rhs_idx].push_back(founder_idx);
							m_idxs_available_rhs.erase(it);
						}
					}
				}
				
				if (edges.end() == edge_range_end)
					break;
				
				edge_range_begin = edge_range_end;
			}
		}
	}
	
	
	void path_mapper::update_string_indices()
	{
		// m_idxs_available_lhs could be used instead.
		std::fill(m_string_idxs_by_founder_lhs.begin(), m_string_idxs_by_founder_lhs.end(), UNASSIGNED_INDEX);
		
		for (auto const &[substring_idx, founder_idx_vec] : rsv::enumerate(m_idxs_by_substring_lhs))
		{
			for (auto const founder_idx : founder_idx_vec)
				m_string_idxs_by_founder_lhs[founder_idx] = substring_idx;
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
	
	
	void sequence_writer::output_subgraph_path(std::size_t const subgraph_idx, substring_index_type const path_idx, std::ostream &os) const
	{
		auto const &ref_positions(m_graph->ref_positions());						// REF positions (0-based) by node number. We use 1-based indexing in order to make summing easier.
		auto const &aln_positions(m_graph->aligned_ref_positions());				// Aligned REF positions by node number. We use 1-based indexing in order to make summing easier.
		auto const &subgraph_start_positions(m_graph->subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		auto const &alt_edge_targets(m_graph->alt_edge_targets());					// ALT edge target nodes by edge number as a concatenated vector.
		auto const &alt_edge_count_csum(m_graph->alt_edge_count_csum());			// Cumulative sum of ALT edges by node number.
		auto const &all_path_edges(m_graph->path_edges());							// Edge numbers (0 for REF edge, 1 for first ALT edge etc.) by path, variant and subgraph number.
		auto const &subgraph_path_edges(all_path_edges[subgraph_idx]);
		// Get a slice that represents the current path. Indices represent the variants (i.e. nodes that have ALT edges).
		auto const &current_path_edges(subgraph_path_edges.column(path_idx));
		
		// Node indices.
		auto const subgraph_begin(subgraph_start_positions[subgraph_idx]);
		auto const subgraph_end(
			1 + subgraph_idx < subgraph_start_positions.size()
			? subgraph_start_positions[1 + subgraph_idx]
			: ref_positions.size() - 2
		);
		
		libbio_assert_eq(os.tellp(), aln_positions[1 + subgraph_begin]);
		
		std::size_t variant_idx(0);
		auto expected_node_idx(subgraph_begin);
		for (auto node_idx(subgraph_begin); node_idx < subgraph_end; ++node_idx)
		{
			libbio_assert_lt(1 + node_idx, alt_edge_count_csum.size());
			auto const alt_edge_start(alt_edge_count_csum[node_idx]);
			auto const alt_edge_limit(alt_edge_count_csum[1 + node_idx]);
			
			if (alt_edge_start == alt_edge_limit)
			{
				// Non-variant node.
				if (node_idx == expected_node_idx)
				{
					++expected_node_idx;
					output_ref(node_idx, expected_node_idx, os);
				}
			}
			else
			{
				// Variant node.
				if (node_idx == expected_node_idx)
				{
					auto const edge(current_path_edges[variant_idx]);
					if (0 == edge)
					{
						++expected_node_idx;
						output_ref(node_idx, expected_node_idx, os);
					}
					else
					{
						auto const alt_idx(alt_edge_start + edge - 1);
						expected_node_idx = alt_edge_targets[alt_idx];
						libbio_assert_lte(expected_node_idx, subgraph_end);
						output_alt(node_idx, expected_node_idx, alt_idx, os);
					}
				}
				++variant_idx;
			}
		}
		
		libbio_assert_eq(os.tellp(), aln_positions[1 + subgraph_end]);
	}
	
	
	void sequence_writer::output_ref_for_subgraph(std::size_t const subgraph_idx, std::ostream &os) const
	{
		auto const &ref_positions(m_graph->ref_positions());
		auto const &subgraph_start_positions(m_graph->subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		
		// Node indices.
		auto const subgraph_begin(subgraph_start_positions[subgraph_idx]);
		auto const subgraph_end(
			1 + subgraph_idx < subgraph_start_positions.size()
			? subgraph_start_positions[1 + subgraph_idx]
			: ref_positions.size() - 2
		);
		
		libbio_assert_eq(os.tellp(), m_graph->aligned_ref_positions()[1 + subgraph_begin]);
		for (auto i(subgraph_begin); i < subgraph_end; ++i)
			output_ref(i, 1 + i, os);
		libbio_assert_eq(os.tellp(), m_graph->aligned_ref_positions()[1 + subgraph_end]);
	}
	
	
	void sequence_writer::output_char_for_subgraph(std::size_t const subgraph_idx, char const cc, std::ostream &os) const
	{
		auto const &subgraph_start_positions(m_graph->subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		auto const &ref_positions(m_graph->ref_positions());
		auto const &aln_positions(m_graph->aligned_ref_positions());
		
		// Node indices.
		auto const subgraph_begin(subgraph_start_positions[subgraph_idx]);
		auto const subgraph_end(
			1 + subgraph_idx < subgraph_start_positions.size()
			? subgraph_start_positions[1 + subgraph_idx]
			: ref_positions.size() - 2
		);
		
		auto const lhs_aln_pos(aln_positions[1 + subgraph_begin]);
		auto const rhs_aln_pos(aln_positions[1 + subgraph_end]);
		libbio_assert_lte(lhs_aln_pos, rhs_aln_pos);
		auto const gap_count(rhs_aln_pos - lhs_aln_pos);
		libbio_assert_eq(os.tellp(), lhs_aln_pos);
		std::fill_n(std::ostream_iterator <char>(os), gap_count, cc);
		libbio_assert_eq(os.tellp(), rhs_aln_pos);
	}
	
	
	void sequence_writer::output_ref(substring_index_type const node_idx, substring_index_type const next_node_idx, std::ostream &os) const
	{
		auto const &ref_positions(m_graph->ref_positions());
		auto const &aln_positions(m_graph->aligned_ref_positions());
		auto const lhs_ref_pos(ref_positions[1 + node_idx]);
		auto const rhs_ref_pos(ref_positions[1 + next_node_idx]);
		auto const lhs_aln_pos(aln_positions[1 + node_idx]);
		auto const rhs_aln_pos(aln_positions[1 + next_node_idx]);
		libbio_assert_lte(lhs_ref_pos, rhs_ref_pos);
		libbio_assert_lte(lhs_aln_pos, rhs_aln_pos);
		auto const ref_len(rhs_ref_pos - lhs_ref_pos);
		auto const aln_len(rhs_aln_pos - lhs_aln_pos);
		libbio_assert_lte(ref_len, aln_len);
		auto const gap_count(aln_len - ref_len);
		auto const ref_sub(m_reference.substr(lhs_ref_pos, ref_len));
		os << ref_sub;
		std::fill_n(std::ostream_iterator <char>(os), gap_count, '-');
	}
	
	
	void sequence_writer::output_alt(substring_index_type const node_idx, substring_index_type const next_node_idx, std::size_t const alt_idx, std::ostream &os) const
	{
		auto const &aln_positions(m_graph->aligned_ref_positions());
		auto const &alt_edge_labels(m_graph->alt_edge_labels());
		auto const lhs_aln_pos(aln_positions[1 + node_idx]);
		auto const rhs_aln_pos(aln_positions[1 + next_node_idx]);
		libbio_assert_lte(lhs_aln_pos, rhs_aln_pos);
		auto const alt_str(alt_edge_labels[alt_idx]);
		auto const aln_len(rhs_aln_pos - lhs_aln_pos);
		libbio_assert_lte(alt_str.size(), aln_len);
		auto const gap_count(aln_len - alt_str.size());
		os << alt_str;
		std::fill_n(std::ostream_iterator <char>(os), gap_count, '-');
	}
}


namespace vcf2multialign {
	
	void founder_sequence_greedy_generator::process_graph_and_output(output_stream_vector &output_files, progress_indicator_delegate &progress_delegate) const
	{
		// Greedy matching is done as follows:
		// 1. Initially, occurring substring numbers in the first segment are assigned to each slot (i.e. output stream).
		//    The remaining slots are marked “free”. (Done with slots_available_lhs.)
		// 2. For each segment pair, the right hand side slots are treated similarly. (Done with slots_available_rhs.)
		// 3. The sorted (by count) path list is then traversed. If a pair of substring numbers is available, they are chosen.
		//    If either lhs or rhs is available and the other side has a free slot, that slot is assigned the corresponding substring index.
		//    If neither lhs or rhs is available but both sides have a free slot, the slots are assigned the corresponding substring indices.
		// 4. Finally, the edges are drawn and the remaining substrings are connected arbitrarily.
		// 5. If a slot is available after all edges have been connected, it is returned to the free list and filled with gap characters.
		
		auto const &sample_paths(m_graph.sample_paths());							// Sample path numbers by sample and subgraph number, vector of vectors.
		auto const &all_path_edges(m_graph.path_edges());
		auto const &subgraph_start_positions(m_graph.subgraph_start_positions());	// Subgraph starting node numbers by subgraph number.
		
		// Before beginning, check whether the first subgraph starts from zero.
		if (subgraph_start_positions.empty())
			return;
		
		sequence_writer sw(m_graph, m_reference);
		
		bool const first_subgraph_starts_from_zero(0 == subgraph_start_positions.front());
		if (!first_subgraph_starts_from_zero)
		{
			// Output REF until the starting position.
			auto const first_subgraph_start_pos(subgraph_start_positions.front());
			for (auto &output_file : output_files)
			{
				for (std::size_t i(0); i < first_subgraph_start_pos; ++i)
					sw.output_ref(i, output_file);
			}
		}
		
		// Count the distinct pairs of paths.
		// A vector is used b.c. we would like to sort by different keys.
		
		path_item_vector path_counts;						// Existing edges in the graph, to be sorted by occurrence count.
		edge_vector edges;									// Edges to be drawn in the founders.
		substring_index_vector substrings_added_to_lhs;		// Substrings added to the lhs segment in the current iteration.
		segment_connector sc(m_founder_count);
		path_mapper pm(m_founder_count);
		
		// Mark the indices in the first segment as assigned for greedy matching.
		std::size_t max_lhs_substring_idx(all_path_edges.front().number_of_columns() - 1);
		libbio_always_assert_lt_msg(max_lhs_substring_idx, m_founder_count, "Given founder count (", m_founder_count, ") is less than the number of distinct substrings in subgraph 1 (", 1 + max_lhs_substring_idx, ").");
		sc.setup(1 + max_lhs_substring_idx);
		pm.setup(1 + max_lhs_substring_idx);
		
		// Handle the segment pairs.
		for (auto const &[lhs_subgraph_idx, pair] : rsv::zip(rsv::ints(0), sample_paths | rsv::sliding(2)))
		{
			// Initialize for the current iteration.
			path_counts.clear();
			edges.clear();
			substrings_added_to_lhs.clear();
			
			// Destructure the pair.
			auto const &lhs_substring_numbers(pair[0]);
			auto const &rhs_substring_numbers(pair[1]);
			libbio_always_assert_eq(lhs_substring_numbers.size(), rhs_substring_numbers.size()); // Verify that the input is valid.
			
			// Iterate over the pairs of paths and count.
			substring_index_type max_rhs_substring_idx(0);
			for (auto const &[lhs_substring_idx, rhs_substring_idx] : rsv::zip(lhs_substring_numbers, rhs_substring_numbers))
			{
				// Try to find an existing path pair.
				auto const res(std::equal_range(path_counts.begin(), path_counts.end(), path_item(lhs_substring_idx, rhs_substring_idx)));
				if (res.first == res.second)
				{
					// Not found (since there were no not-less-than elements). Append to the end and rotate the greater-than range.
					auto const first_greater_than_idx(std::distance(path_counts.begin(), res.second)); // Index of the first greater-than element.
					path_counts.emplace_back(lhs_substring_idx, rhs_substring_idx, 1); // Invalidates iterators.
					
					auto const begin(path_counts.begin() + first_greater_than_idx);
					auto const end(path_counts.end());
					std::rotate(begin, end - 1, end); // There is at least one element b.c. we just emplace_back’d one.
				}
				else
				{
					// Found.
					libbio_assert_eq(1, std::distance(res.first, res.second)); // Check that there is in fact only one matching item.
					res.first->count += 1;
				}
				
				max_rhs_substring_idx = lb::max_ct(max_rhs_substring_idx, rhs_substring_idx);
			}
			
			libbio_always_assert_lt_msg(max_rhs_substring_idx, m_founder_count, "Given founder count (", m_founder_count, ") is less than the number of distinct substrings in subgraph ", 1 + lhs_subgraph_idx, " (", 1 + max_rhs_substring_idx, ").");
			
			// Sort by count and add edges in descending count order.
			std::sort(path_counts.begin(), path_counts.end(), [](auto const &lhs, auto const &rhs) -> bool {
				return lhs.count < rhs.count;
			});
			sc.make_edges(path_counts, 1 + max_rhs_substring_idx, edges, substrings_added_to_lhs);
			
			// Associate the edges with the founders i.e. output streams.
			// assign_edges_to_founders requires the edges to be sorted by lhs_idx.
			std::sort(edges.begin(), edges.end());
			pm.add_substrings(substrings_added_to_lhs);
			pm.assign_edges_to_founders(edges);
			pm.update_string_indices();
			
			// Output the founders.
			for (auto const &[founder_idx, substring_idx] : rsv::enumerate(pm.string_indices_by_founder()))
			{
				auto &os(output_files[founder_idx]);
				if (UNASSIGNED_INDEX == substring_idx)
					sw.output_n_for_subgraph(lhs_subgraph_idx, os);
				else
					sw.output_subgraph_path(lhs_subgraph_idx, substring_idx, os);
			}
			
			if (m_output_reference)
			{
				auto &os(output_files.back());
				sw.output_ref_for_subgraph(lhs_subgraph_idx, os);
			}
			
			pm.end_subgraph();
			progress_delegate.advance();
			max_lhs_substring_idx = max_rhs_substring_idx;
		}
		
		// Handle the last subgraph and output the founders.
		{
			pm.update_string_indices();
			auto const last_subgraph_idx(m_graph.subgraph_count() - 1);
			for (auto const &[founder_idx, substring_idx] : rsv::enumerate(pm.string_indices_by_founder()))
			{
				auto &os(output_files[founder_idx]);
				if (UNASSIGNED_INDEX == substring_idx)
					sw.output_n_for_subgraph(last_subgraph_idx, os);
				else
					sw.output_subgraph_path(last_subgraph_idx, substring_idx, os);
				os << std::flush;
			}
			
			if (m_output_reference)
			{
				auto &os(output_files.back());
				sw.output_ref_for_subgraph(last_subgraph_idx, os);
				os << std::flush;
			}
		}
		progress_delegate.advance();
	}
	
	
	void founder_sequence_greedy_generator::output_sequences()
	{
		try
		{
			// Setup the progress indicator.
			this->install_progress_indicator();
			progress_indicator_delegate progress_delegate(m_graph.subgraph_count());
			this->progress_indicator().log_with_progress_bar("\t", progress_delegate);
			
			// Create a string view from the reference.
			std::string_view const reference_sv(m_reference.data(), m_reference.size());
			
			auto const output_count(m_founder_count + m_output_reference);
			dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << m_founder_count << " sequences will be written.\n"; });
			
			// Don’t use chunks for now b.c. that would complicate things too much and not writing everything simultaneously
			// is more useful with predicted sequences (not founders).
			
			output_stream_vector output_files(output_count);
			
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(m_may_overwrite ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			
			for (std::size_t i(0); i < m_founder_count; ++i)
				open_founder_output_file(i, output_files[i], mode);
			
			if (m_output_reference)
				lb::open_file_for_writing("REF", output_files.back(), mode);
			
			// Generate the sequences.
			process_graph_and_output(output_files, progress_delegate);
			
			dispatch_async(dispatch_get_main_queue(), ^{
				lb::log_time(std::cerr);
				std::cerr << "Done.\n"; // FIXME: log statistics?
				this->finish_mt();
			});
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			this->log_assertion_failure_and_exit(exc);
		}
		catch (std::exception const &exc)
		{
			this->log_exception_and_exit(exc);
		}
		catch (...)
		{
			this->log_unknown_exception_and_exit();
		}
	}
}
