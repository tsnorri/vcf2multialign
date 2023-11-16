/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>							// std::lower_bound, std::max, std::reverse, std::upper_bound
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/subrange.hpp>
#include <vcf2multialign/find_cut_positions.hh>
#include <vcf2multialign/pbwt.hh>

namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef v2m::pbwt_context <
		v2m::variant_graph::sample_type,
		v2m::variant_graph::edge_type,
		v2m::variant_graph::ploidy_type
	>										pbwt_context_type;
	
	
	// We calculate positions by edge numbers due to the fact that a path using a given edge is a binary property
	// and hence we would like to maintain the pBWT divergence values for them.
	struct cut_position
	{
		typedef v2m::variant_graph::node_type	node_type;
		typedef v2m::variant_graph::edge_type	edge_type;
		typedef v2m::variant_graph::ploidy_type	count_type;
		
		edge_type	edge{};										// The first edge in the node by which to cut.
		edge_type	prev_edge{v2m::variant_graph::EDGE_MAX};
		node_type	node{};
		count_type	score{};
		
		void update_if_needed(count_type const eq_class_count, cut_position const &prev_cut);
	};
	
	typedef std::vector <cut_position> cut_position_vector;
	
	struct cut_position_cmp
	{
		typedef v2m::variant_graph::edge_type	edge_type;
		
		bool operator()(cut_position const &lhs, cut_position const &rhs) const { return lhs.edge < rhs.edge; }
		bool operator()(cut_position const &lhs, edge_type const rhs) const { return lhs.edge < rhs; }
		bool operator()(edge_type const lhs, cut_position const &rhs) const { return lhs < rhs.edge; }
	};
	
	void cut_position::update_if_needed(count_type const eq_class_count, cut_position const &prev_cut)
	{
		auto const candidate_score{std::max(eq_class_count, prev_cut.score)};
		if (candidate_score < score)
		{
			score = candidate_score;
			prev_edge = prev_cut.edge;
		}
	}
}


namespace vcf2multialign {
	
	// Find cut positions in the graph minimising the block height.
	// The algorithm uses pBWT to determine the number of equivalence classes
	// of the sequence segments between candidate cut positions. To use the
	// binary alphabet version of pBWT, we consider each ALT edge separately
	// instead of each node. A node is a candidate cut position if it is an
	// endpoint of a bridge.
	//
	// The algorithm works as follows. In addition to the a and d arrays of the
	// pBWT, we maintain a map of divergence value counts.
	//	– When we arrive at a node, we check if it is a candidate cut position.
	//		– If this is the case, we calculate the scores of the subgraphs ending at
	//		  said position and pick the best one.
	//		– This is done by iterating over the (at most m) divergence values, picking
	//		  the leftmost unhandled cut position the (edge) index of which is greater than
	//		  the one that corresponds to the divergence value and calculating the score.
	//		– The divergence values are handled from right to left, i.e. that the
	//		  smallest number of equivalence classes is considered first. Each candidate
	//		  cut position needs to be considered at most once, since the score of the
	//		  graph segment being calculated will increase when the number of equivalence
	//		  classes is increased.
	//		– Finally, we consider the case where the current subgraph extends beyond the
	//		  leftmost divergence value. (This is particularly helpful when the aligned length
	//		  of the current subgraph is less than min_length.)
	//	– Before leaving the node, we update the pBWT values for each ALT edge separately.
	cut_position_score_type find_initial_cut_positions_lambda_min(
		variant_graph const &graph,
		variant_graph::edge_type const min_distance,
		std::vector <variant_graph::position_type> &out_cut_positions
	)
	{
		out_cut_positions.clear();
	
		variant_graph::node_type rightmost_seen_alt_edge_target{};
		variant_graph_walker walker(graph);
		variant_graph::edge_type edge_idx{};
	
		pbwt_context_type pbwt_ctx(graph.path_count());
		cut_position_vector cut_positions;
	
		// Divergence value counts without the “k + 1” count.
		auto const divergence_value_counts_reversed([&pbwt_ctx]{
			auto const &dvc(pbwt_ctx.divergence_value_counts);
			auto const it(dvc.begin());
			auto end(dvc.end());
			--end;
			return rsv::reverse(ranges::subrange(it, end));
		});
	
		while (walker.advance())
		{
			// Check if the current node is a potential cut position.
			if (rightmost_seen_alt_edge_target <= walker.node())
			{
				// Process the divergence values from the smallest number of equivalence classes.
				// Since we use edge_idx + walker.alt_edge_count() here, we need to use std::upper_bound()
				// to find the correct solution later.
				auto const cut_pos_id(edge_idx + walker.alt_edge_count());
				auto &current_cut(cut_positions.emplace_back(cut_pos_id, variant_graph::EDGE_MAX, walker.node(), graph.path_count()));
			
				auto const cut_pos_begin(cut_positions.begin());
				auto const cut_pos_end(cut_positions.end());
				auto cut_pos_rb(cut_pos_end);
				auto eq_class_count(pbwt_ctx.divergence_value_counts.rbegin()->second);
				for (auto const &[div_edge_idx, div_count] : divergence_value_counts_reversed())
				{
					// Find the leftmost cut position greater than div_edge_idx.
					// We only need to check said position once b.c. the number of equivalence classes,
					// as determined from the divergence values, incereases as we iterate over said values.
					auto const it(std::upper_bound(cut_pos_begin, cut_pos_rb, div_edge_idx, cut_position_cmp{}));
					if (it != cut_pos_end)
					{
						cut_pos_rb = it;
					
						// Check if the distance to the previous node is at least min_distance.
						// FIXME: alternatively we could use the minimum path length between said nodes. Calculating it is more difficult, though.
						if (min_distance <= graph.aligned_length(it->node, walker.node()))
						{
							// Check if we can improve the score.
							current_cut.update_if_needed(eq_class_count, *it);
						}
					}
				
					eq_class_count += div_count;
				}
			
				// Check the cut position immediately to the left from cut_pos_rb.
				if (cut_pos_begin != cut_pos_rb)
				{
					--cut_pos_rb;
					current_cut.update_if_needed(eq_class_count, *cut_pos_rb);
				}
			}
		
			// Handle the edges.
			for (auto const dst_node : walker.alt_edge_targets())
			{
				pbwt_ctx.swap_vectors();
				pbwt_ctx.update_divergence(graph.paths_by_edge_and_chrom_copy.column(edge_idx), edge_idx);
				++edge_idx;
				rightmost_seen_alt_edge_target = std::max(rightmost_seen_alt_edge_target, dst_node);
			}
		}
	
		// Copy the solution if possible.
		if (cut_positions.empty())
			return CUT_POSITION_SCORE_MAX;
		
		
		{
			cut_position_vector::const_iterator it(cut_positions.end() - 1);
			auto const retval(it->score);
			while (true)
			{
				auto const node(it->node);
				out_cut_positions.push_back(node);
				auto const prev_edge(it->prev_edge);
				if (variant_graph::EDGE_MAX == prev_edge)
					break;
			
				it = std::lower_bound(cut_positions.begin(), cut_positions.end(), prev_edge, cut_position_cmp{});
			}
		
			std::reverse(out_cut_positions.begin(), out_cut_positions.end());
			
			// Add a sentinel to make processing easier.
			out_cut_positions.push_back(variant_graph::NODE_MAX);
			
			return retval;
		}
	}
}
