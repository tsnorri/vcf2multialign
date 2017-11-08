/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_MERGE_SUBGRAPH_PATHS_TASK_HH
#define VCF2MULTIALIGN_TASKS_MERGE_SUBGRAPH_PATHS_TASK_HH

#include <lemon/full_graph.h>
#include <lemon/matching.h>
#include <mutex>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/reduced_subgraph.hh>
#include <vcf2multialign/status_logger.hh>
#include <vcf2multialign/tasks/task.hh>


namespace vcf2multialign {
	
	struct merge_subgraph_paths_task_delegate;
	
	
	class merge_subgraph_paths_task : public task
	{
	public:
		typedef int32_t																weight_type;
		
	protected:
		typedef lemon::FullBpGraph													graph_type;
		typedef int32_t																edge_cost_type;
		typedef std::vector <edge_cost_type>										edge_cost_vector_type;
		typedef graph_type::EdgeMap <edge_cost_type>								edge_cost_map_type;
		typedef lemon::MaxWeightedPerfectMatching <graph_type, edge_cost_map_type>	matching_type;	// No bipartite algorithms in Lemon 1.3.1.
		
	protected:
		merge_subgraph_paths_task_delegate	*m_delegate{nullptr};
		status_logger						*m_status_logger{nullptr};
		reduced_subgraph const				*m_left_subgraph{nullptr};
		reduced_subgraph const				*m_right_subgraph{nullptr};
		std::size_t							m_lhs_idx{0};
		std::size_t							m_path_count{0};
		
	public:
		merge_subgraph_paths_task() = default;
		
		merge_subgraph_paths_task(
			merge_subgraph_paths_task_delegate &delegate,
			class status_logger &status_logger,
			reduced_subgraph const &left,
			reduced_subgraph const &right,
			std::size_t const lhs_idx,
			std::size_t const path_count
		):
			m_delegate(&delegate),
			m_status_logger(&status_logger),
			m_left_subgraph(&left),
			m_right_subgraph(&right),
			m_lhs_idx(lhs_idx),
			m_path_count(path_count)
		{
		}
		
		std::size_t left_subgraph_index() const { return m_lhs_idx; }
		
		virtual void execute() override;
		
	protected:
		void calculate_edge_weight(
			graph_type const &graph,
			edge_cost_vector_type &edge_costs,
			std::size_t const li,
			std::size_t const ri
		);
			
		weight_type find_minimum_cost_matching(
			graph_type const &graph,
			edge_cost_map_type const &edge_costs,
			std::vector <reduced_subgraph::path_index> &matchings
		);
	};
	
	
	struct merge_subgraph_paths_task_delegate
	{
		virtual ~merge_subgraph_paths_task_delegate() {}
		virtual void task_did_calculate_edge_weight(merge_subgraph_paths_task &task) {};
		virtual void task_did_finish(
			merge_subgraph_paths_task &task,
			std::vector <reduced_subgraph::path_index> &&matchings,
			merge_subgraph_paths_task::weight_type const matching_weight
		) = 0;
	};
}


#endif
