/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/tasks/merge_subgraph_paths_task.hh>


namespace vcf2multialign {
	
	void merge_subgraph_paths_task::calculate_edge_weight(
		graph_type const &graph,
		edge_cost_vector_type &edge_costs,
		std::size_t const li,
		std::size_t const ri
	)
	{
		// Try to make sure that the weight is within data type limits.
		auto const inverse_weight(edge_weight(*m_left_subgraph, *m_right_subgraph, li, ri));
		weight_type const weight(-inverse_weight);
		always_assert(-weight == inverse_weight);
		
		// Store the weight.
		auto const lhs(graph.redNode(li));
		auto const rhs(graph.blueNode(ri));
		auto const edge(graph.edge(lhs, rhs));
		auto const edge_id(graph_type::id(edge));
		always_assert(edge_id < edge_costs.size());
		edge_costs[edge_id] = weight;
		
		m_delegate->task_did_calculate_edge_weight(*this);
	}
	
	
	auto merge_subgraph_paths_task::find_minimum_cost_matching(
		graph_type const &graph,
		edge_cost_map_type const &edge_costs,
		std::vector <reduced_subgraph::path_index> &matchings
	) -> weight_type
	{
		matching_type matching(graph, edge_costs);
		matching.run();

		// Iterate red nodes in the graph, find their mates and store the pairs in target_paths.
		auto const matching_weight(matching.matchingWeight());
		for (std::size_t i(0), count(graph.redNum()); i < count; ++i)
		{
			auto const lhs(graph.redNode(i));
			always_assert(lhs != lemon::INVALID);
			auto const node(matching.mate(lhs));
			always_assert(node != lemon::INVALID);
			auto const rhs(graph.asBlueNode(node));
			always_assert(rhs != lemon::INVALID);
			
			auto const li(graph.index(lhs));
			auto const ri(graph.index(rhs));
			always_assert(li < matchings.size());
			always_assert(ri < matchings.size());
			matchings[li] = ri;
		}
		
		return matching_weight;
	}
	
	
	void merge_subgraph_paths_task::execute()
	{
		always_assert(m_left_subgraph);
		always_assert(m_right_subgraph);
		always_assert(m_path_count);
		
		dispatch_queue_t queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
		dispatch_ptr <dispatch_group_t> group(dispatch_group_create());
		
		// Use a std::vector since contents of its contained objects in different elements
		// may be modified concurrently without data races.
		graph_type graph(m_path_count, m_path_count);
		edge_cost_map_type edge_cost_map(graph);
		
		{
			edge_cost_vector_type edge_costs(1 + graph.maxEdgeId(), 0);
			for (std::size_t i(0); i < m_path_count; ++i)
			{
				for (std::size_t j(0); j < m_path_count; ++j)
				{
					dispatch_group_async_fn(*group, queue, [this, &graph, &edge_costs, i, j](){
						calculate_edge_weight(graph, edge_costs, i, j);
					});
				}
			}
			
			dispatch_group_wait(*group, DISPATCH_TIME_FOREVER);
			
			// Copy the values to edge_cost_map.
			std::size_t i(0);
			for (auto const weight : edge_costs)
			{
				auto const edge(graph_type::edgeFromId(i));
				edge_cost_map.set(edge, weight);
				++i;
			}
		}
		
		std::vector <reduced_subgraph::path_index> matchings(
			m_path_count,
			std::numeric_limits <reduced_subgraph::path_index>::max()
		);

		auto const matching_weight(find_minimum_cost_matching(graph, edge_cost_map, matchings));
		m_delegate->task_did_finish(*this, std::move(matchings), matching_weight);
	}
}
