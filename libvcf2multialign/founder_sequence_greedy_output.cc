/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>							// std::fill, std::sort
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <libbio/assert.hh>
#include <libbio/file_handling.hh>
#include <libbio/int_vector.hh>					// lb::bit_vector
#include <libbio/matrix.hh>						// lb::matrix
#include <map>									// std::multimap
#include <range/v3/view/all.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/zip.hpp>
#include <vcf2multialign/find_cut_positions.hh>
#include <vcf2multialign/pbwt.hh>
#include <vcf2multialign/output.hh>
#include <vcf2multialign/sequence_writer.hh>
#include <utility>								// std::swap
#include <vector>

namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef v2m::founder_sequence_greedy_output::ploidy_type	ploidy_type;
	
	typedef v2m::pbwt_context <
		v2m::variant_graph::sample_type,
		v2m::variant_graph::edge_type,
		ploidy_type
	>															pbwt_context_type;
	
	
	struct joined_path_eq_class
	{
		ploidy_type	lhs_rep{};
		ploidy_type	rhs_rep{};
		ploidy_type	size{};
		
		joined_path_eq_class(ploidy_type const lhs_rep_, ploidy_type const rhs_rep_):
			lhs_rep(lhs_rep_),
			rhs_rep(rhs_rep_)
		{
		}
		
		bool operator<(joined_path_eq_class const &other) const { return size < other.size; }
	};
	
	
	struct reference_sequence_writing_delegate final : public v2m::sequence_writing_delegate
	{
		void handle_node(variant_graph const &graph, node_type const node) override {}
	};
	
	
	class founder_sequence_writing_delegate final : public v2m::sequence_writing_delegate
	{
	public:
		typedef v2m::variant_graph::position_type	position_type;
		typedef v2m::founder_sequence_greedy_output	output_type;
		typedef output_type::ploidy_matrix			ploidy_matrix;
		typedef ploidy_matrix::const_slice_type		ploidy_matrix_const_slice;
		typedef output_type::cut_position_vector	cut_position_vector;
		
	private:
		ploidy_matrix_const_slice const	m_assigned_samples;
		cut_position_vector const		&m_cut_positions;
		position_type					m_cut_pos_index{1};
		
	public:
		founder_sequence_writing_delegate(
			ploidy_matrix_const_slice &&assigned_samples,
			cut_position_vector const &cut_positions
		):
			v2m::sequence_writing_delegate(),
			m_assigned_samples(std::move(assigned_samples)),
			m_cut_positions(cut_positions)
		{
			libbio_assert(!m_cut_positions.empty());
			libbio_assert_eq(0, m_cut_positions.front());
		}
		
		
		void handle_node(variant_graph const &graph, node_type const node) override
		{
			libbio_assert_lte(node, m_cut_positions[m_cut_pos_index]);
			if (node == m_cut_positions[m_cut_pos_index])
			{
				chromosome_copy_index = m_assigned_samples[m_cut_pos_index - 1];
				++m_cut_pos_index;
			}
		}
	};
}


namespace vcf2multialign {
	
	void founder_sequence_greedy_output::load_cut_positions(char const *path)
	{
		lb::file_istream is;
		lb::open_file_for_reading(path, is);
		cereal::PortableBinaryInputArchive archive(is);
		archive(m_cut_positions);
	}
	
	
	void founder_sequence_greedy_output::output_cut_positions(char const *path)
	{
		lb::file_ostream os;
		lb::open_file_for_writing(path, os, lb::writing_open_mode::CREATE);
		cereal::PortableBinaryOutputArchive archive(os);
		archive(m_cut_positions);
	}
	
	
	bool founder_sequence_greedy_output::find_cut_positions(
		variant_graph const &graph,
		variant_graph::position_type const min_dist
	)
	{
		auto const score(find_initial_cut_positions_lambda_min(graph, min_dist, m_cut_positions.cut_positions));
		if (CUT_POSITION_SCORE_MAX == score)
			return false;
		
		m_cut_positions.min_distance = min_dist;
		m_cut_positions.score = score;
		return true;
	}
	
	
	bool founder_sequence_greedy_output::find_matching(variant_graph const &graph, ploidy_type const founder_count)
	{
		// We re-calculate the pBWT in order to determine the equivalence class representatives
		// of paths between adjacent cut positions. When we have a pair of such blocks (and lists
		// of representatives), we re-use the just calculated pBWT to determine the equivalence
		// classes of paths from the left cutting position of the pair to the right one.
		// Finally we use the sizes of the resulting equivalence classes in the mathcing.
		
		if (m_cut_positions.cut_positions.size() < 2)
			return false;
		
		libbio_assert_eq(0, m_cut_positions.cut_positions.front());
		libbio_assert_eq(graph.node_count() - 1, m_cut_positions.cut_positions.back());
		
		m_assigned_samples.clear();
		m_assigned_samples.resize(m_cut_positions.cut_positions.size(), founder_count); // Founders in columns.
		std::fill(m_assigned_samples.begin(), m_assigned_samples.end(), PLOIDY_MAX);
		
		std::multimap <ploidy_type, ploidy_type> assignments_by_eq_class;
		lb::bit_vector reserved_assignments(graph.total_chromosome_copies(), 0);
		std::vector <ploidy_type> arbitrarily_connected_rhs;
		
		variant_graph_walker walker(graph);
		variant_graph::edge_type edge_idx{};
		variant_graph::edge_type prev_cut_edge_idx{};
		variant_graph::edge_type cut_pair_edge_idx{};
		
		std::vector <ploidy_type> lhs_eq_classes(graph.total_chromosome_copies(), PLOIDY_MAX);
		std::vector <ploidy_type> rhs_eq_classes(graph.total_chromosome_copies(), PLOIDY_MAX);
		ploidy_type lhs_distinct_eq_classes{};
		ploidy_type rhs_distinct_eq_classes{};
		std::vector <joined_path_eq_class> joined_path_eq_classes;
		
		// m_cut_positions has a value for the sink node.
		auto cut_pos_it(m_cut_positions.cut_positions.begin());
		++cut_pos_it; // Node zero.
		
		pbwt_context_type pbwt_ctx(graph.total_chromosome_copies());
		
		// Handle the rest.
		for (variant_graph::position_type cut_pos_idx{}; walker.advance(); ++cut_pos_idx)
		{
			libbio_assert_neq(cut_pos_it, m_cut_positions.cut_positions.end());
			
			auto const node(walker.node());
			libbio_assert_lte(node, *cut_pos_it);
			
			// Check if we are at a cut position.
			if (node == *cut_pos_it)
			{
				{
					using std::swap;
					swap(lhs_eq_classes, rhs_eq_classes);
					swap(lhs_distinct_eq_classes, rhs_distinct_eq_classes);
				}
				
				{
					// Determine the rhs. and the joined eq. classes.
					// Note that due to how these are determined, the class representatives
					// are not interchangeable between blocks (separated by cut positions).
					ploidy_type rep{PLOIDY_MAX};
					rhs_distinct_eq_classes = 0;
					for (auto const [aa, dd] : rsv::zip(pbwt_ctx.permutation, pbwt_ctx.divergence))
					{
						// Check if the current entry begins a new equivalence class.
						if (prev_cut_edge_idx < dd)
						{
							rep = aa;
							++rhs_distinct_eq_classes;
						}
						
						// Store for the next cut position.
						rhs_eq_classes[aa] = rep;
						
						// We rely on the branch predictor to take care of this.
						if (1 < cut_pos_idx)
						{
							if (cut_pair_edge_idx < dd)
								joined_path_eq_classes.emplace_back(lhs_eq_classes[aa], rep);
						
							libbio_assert(!joined_path_eq_classes.empty());
							++joined_path_eq_classes.back().size;
						}
					}
				}
				
				if (2 <= cut_pos_idx)
				{
					// Sort by the size. (The smallest will be the first.)
					std::sort(joined_path_eq_classes.begin(), joined_path_eq_classes.end());
					
					if (2 == cut_pos_idx)
					{
						// Second cut position; initial assignment.
						
						auto remaining_founders(founder_count);
						remaining_founders -= lhs_distinct_eq_classes;
						
						ploidy_type founder_idx{};
						
						auto const do_assign([&](joined_path_eq_class const &eq_class){
							assignments_by_eq_class.emplace(eq_class.lhs_rep, founder_idx);
							m_assigned_samples(cut_pos_idx - 1, founder_idx) = eq_class.lhs_rep;
							++founder_idx;
						});
						
						for (auto const &eq_class : rsv::reverse(joined_path_eq_classes))
						{
							auto ref(reserved_assignments[eq_class.lhs_rep]);
							static_assert(ref.is_reference()); // Sanity check.
							if (ref)
							{
								// Already seen the eq. class (i.e. used the reserved slot);
								// Try to get a new one.
								if (remaining_founders)
								{
									--remaining_founders;
									do_assign(eq_class);
								}
							}
							else
							{
								// Mark seen and place a copy.
								ref |= 0x1;
								do_assign(eq_class);
							}
						}
						
						// We would like to have the invariant that all founders have an
						// assigned eq. class.
						while (true)
						{
							for (auto const &eq_class : rsv::reverse(joined_path_eq_classes))
							{
								if (!remaining_founders)
									goto handle_subsequent_assignment;
							
								--remaining_founders;
								do_assign(eq_class);
							}
						}
					}
					
					// Handle the subsequent assignment as follows.
					// 1.	Sp. the eq. class on the right is a reserved one. We check if there is a suitable assignment
					//		on the left. If this is true, we assign. If not, we do nothing.
					// 2.	Sp. not but there is a founder available on the right hand side. Again we check if there
					//		is a suitable assignment on the left. If this is true, we assign. If not, we do nothing.
					// 3.	We continue (2) until all the founders have been used or no assignments were made.
					// 4.	We go through the distinct eq. classes on the right hand side and connect arbitrarily.
					// 5.	We assign rhs sequences to the remaining founders and connect arbitrarily.
					
				handle_subsequent_assignment:
					{
						std::fill(reserved_assignments.word_begin(), reserved_assignments.word_end(), 0);
						arbitrarily_connected_rhs.clear();
						
						auto remaining_founders(founder_count);
						remaining_founders -= rhs_distinct_eq_classes;
						
						auto const try_assign([&](joined_path_eq_class const &eq_class) -> bool {
							auto const it(assignments_by_eq_class.find(eq_class.lhs_rep));
							if (assignments_by_eq_class.end() != it)
							{
								// Found a suitable assignment.
								auto const founder_idx(it->second);
								assignments_by_eq_class.erase(it);
								m_assigned_samples(cut_pos_idx, founder_idx) = eq_class.rhs_rep;
								return true;
							}
							
							return false;
						});
						
						auto const assign_arbitrary([&](ploidy_type const rhs_rep){
							libbio_assert(!assignments_by_eq_class.empty());
							auto const it(assignments_by_eq_class.begin());
							auto const founder_idx(it->second);
							assignments_by_eq_class.erase(it);
							m_assigned_samples(cut_pos_idx, founder_idx) = rhs_rep;
						});
						
						// 1, 2, 3.
						{
							bool is_first{true};
							bool did_assign{false};
							while (true)
							{
								for (auto const &eq_class : rsv::reverse(joined_path_eq_classes))
								{
									auto ref(reserved_assignments[eq_class.rhs_rep]);
									static_assert(ref.is_reference()); // Sanity check.
									if (ref)
									{
										// Already seen.
										if (remaining_founders)
										{
											if (try_assign(eq_class))
											{
												did_assign = true;
												--remaining_founders;
											}
										}
										else if (!is_first) // Small optimisation.
										{
											goto continue_subsequent_assignment;
										}
									}
									else
									{
										// Mark seen and place a copy.
										// Not adding to arbitrarily_connected_rhs on success is actually a small
										// optimisation since we re-check reserved_assignments in (4) anyway.
										if (try_assign(eq_class))
											ref |= 0x1;
										else
											arbitrarily_connected_rhs.push_back(eq_class.rhs_rep);
									}
								}
								
								if (!remaining_founders)
									break;
								
								if (is_first)
								{
									is_first = false;
									continue;
								}
								
								if (!did_assign)
									break;
							}
						}
						
					continue_subsequent_assignment:
						// 4.
						for (auto const rhs_rep : arbitrarily_connected_rhs)
						{
							auto ref(reserved_assignments[rhs_rep]);
							static_assert(ref.is_reference()); // Sanity check.
							if (!ref)
							{
								assign_arbitrary(rhs_rep);
								ref |= 0x1;
							}
						}
						
						// 5.
						while (!assignments_by_eq_class.empty())
						{
							for (auto const &eq_class : rsv::reverse(joined_path_eq_classes))
							{
								if (assignments_by_eq_class.empty())
									goto finish_assignment;
								
								assign_arbitrary(eq_class.rhs_rep);
							}
						}
						
						// Update assignments_by_eq_class to reflect the new state.
					finish_assignment:
						{
							assignments_by_eq_class.clear();
							auto const row(m_assigned_samples.const_row(cut_pos_idx));
							for (auto const &[idx, eq_class] : row | rsv::enumerate)
								assignments_by_eq_class.insert({eq_class, idx});
						}
					}
				}
				
				++cut_pos_it;
				cut_pair_edge_idx = prev_cut_edge_idx;
				prev_cut_edge_idx = edge_idx;
			}
			
			// Handle the edges.
			for (auto const dst_node : walker.alt_edge_targets())
			{
				pbwt_ctx.swap_vectors();
				pbwt_ctx.update_divergence(graph.paths_by_edge_and_chrom_copy.column(edge_idx), edge_idx);
				++edge_idx;
			}
		}
		
		return true;
	}
	
	
	void founder_sequence_greedy_output::output_sequences_a2m(sequence_type const &ref_seq, variant_graph const &graph, lb::file_ostream &stream)
	{
		typedef variant_graph::ploidy_type	ploidy_type;
		
		{
			stream << ">REF\n";
			reference_sequence_writing_delegate delegate;
			output_sequence(ref_seq, graph, stream, delegate);
			stream << '\n';
			m_delegate->handled_sequences(1);
		}
		
		ploidy_type const col_count(m_assigned_samples.number_of_columns());
		for (auto const col_idx : rsv::iota(ploidy_type(0), col_count))
		{
			m_delegate->will_handle_founder_sequence(col_idx);
			
			stream << '>' << (1 + col_idx) << '\n';
			founder_sequence_writing_delegate delegate(m_assigned_samples.const_column(col_idx), m_cut_positions.cut_positions);
			output_sequence(ref_seq, graph, stream, delegate);
			stream << '\n';
			
			m_delegate->handled_sequences(2 + col_idx);
		}
	}
	
	
	void founder_sequence_greedy_output::output_separate(sequence_type const &ref_seq, variant_graph const &graph)
	{
		typedef variant_graph::ploidy_type	ploidy_type;
		
		{
			reference_sequence_writing_delegate delegate;
			output_sequence_file(ref_seq, graph, "REF", delegate);
		}
		
		ploidy_type const col_count(m_assigned_samples.number_of_columns());
		for (auto const col_idx : rsv::iota(ploidy_type(0), col_count))
		{
			m_delegate->will_handle_founder_sequence(col_idx);
			
			auto const dst_name{std::to_string(1 + col_idx)};
			founder_sequence_writing_delegate delegate(m_assigned_samples.const_column(col_idx), m_cut_positions.cut_positions);
			output_sequence_file(ref_seq, graph, dst_name.data(), delegate);
		}
	}
}
