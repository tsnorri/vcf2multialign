/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_REDUCE_SAMPLES_TASK_HH
#define VCF2MULTIALIGN_TASKS_REDUCE_SAMPLES_TASK_HH

#include <numeric>
#include <vcf2multialign/reduced_subgraph.hh>
#include <vcf2multialign/tasks/merge_subgraph_paths_task.hh>
#include <vcf2multialign/tasks/read_subgraph_variants_task.hh>
#include <vcf2multialign/tasks/sequence_writer_task.hh>
#include <vcf2multialign/tasks/task.hh>


namespace vcf2multialign { namespace detail {
	class reduce_samples_progress_counter : public status_logger_delegate
	{
	protected:
		std::atomic_size_t		m_current_step{0};
		std::size_t				m_step_count{0};
		
		// For debugging.
		std::atomic_size_t		m_rsv_steps{0};
		std::atomic_size_t		m_edge_weight_steps{0};
		std::atomic_size_t		m_merge_tasks{0};
		
	public:
		void calculate_step_count(
			std::size_t const valid_record_count,
			std::size_t const path_count,
			std::size_t const merge_task_count
		);
		
		inline void reset_step_count(std::size_t const step_count) { m_current_step = 0; m_step_count = step_count; }
		
		inline void read_subgraph_variants_task_did_handle_variant() { ++m_current_step; --m_rsv_steps; }
		inline void merge_subgraph_paths_task_did_finish() { ++m_current_step; --m_merge_tasks; }
		inline void merge_subgraph_paths_task_did_calculate_edge_weight() { ++m_current_step; --m_edge_weight_steps; }
		inline void sequence_writer_task_did_handle_variant() { ++m_current_step; }
		
		// status_logger_delegate
		virtual std::size_t step_count() const override { return m_step_count; }
		virtual std::size_t current_step() const override { return m_current_step; }
	};
}}


namespace vcf2multialign {
	
	class reduce_samples_task;
	
	
	struct reduce_samples_task_delegate
	{
		virtual void store_and_execute(std::unique_ptr <task> &&task) = 0;
		virtual void task_did_finish(reduce_samples_task &task) = 0;
	};
	
	
	class reduce_samples_task :
		public task,
		public merge_subgraph_paths_task_delegate,
		public read_subgraph_variants_task_delegate,
		public sequence_writer_task_delegate
	{
	public:
		typedef haplotype_map <channel_ostream>			haplotype_map_type;
		typedef std::vector <reduced_subgraph>			subgraph_vector;
		typedef std::vector <
			std::vector <
				reduced_subgraph::path_index
			>
		>												path_matchings_type;
		typedef merge_subgraph_paths_task::weight_type	weight_type;
		typedef std::atomic <weight_type>				atomic_weight_type;
		
		
	protected:
		dispatch_ptr <dispatch_semaphore_t>			m_subgraph_semaphore;
		dispatch_ptr <dispatch_semaphore_t>			m_write_semaphore;

		std::mutex									m_subgraph_mutex;
		subgraph_vector								m_subgraphs;
		std::vector <bool>							m_subgraph_bitmap;
		path_matchings_type							m_path_matchings;
		std::vector <std::size_t>					m_path_permutation;
		haplotype_map_type							m_haplotypes;
		subgraph_vector::const_iterator				m_subgraph_iterator{};
		vcf_reader									m_reader;
		detail::reduce_samples_progress_counter		m_progress_counter;
		
		std::atomic_int32_t							m_running_subgraph_tasks{0};
		std::atomic_int32_t							m_running_merge_tasks{0};
		
		reduce_samples_task_delegate				*m_delegate{nullptr};
		generate_configuration const				*m_generate_config{nullptr};
		struct logger								*m_logger{nullptr};
		preprocessing_result const					*m_preprocessing_result{nullptr};
		mmap_handle const							*m_vcf_input_handle{nullptr};
		std::atomic_size_t							m_remaining_merge_tasks{0};
		atomic_weight_type							m_matching_weight{0};
		std::size_t									m_record_count{0};
		std::size_t									m_sample_ploidy_sum{0};
		
	public:
		reduce_samples_task() = default;
		
		reduce_samples_task(
			reduce_samples_task_delegate &delegate,
			generate_configuration const &generate_config,
			struct logger &logger,
			vcf_reader &&reader,
			mmap_handle const &vcf_input_handle,
			preprocessing_result const &result,
			std::size_t const hw_concurrency,
			std::size_t const record_count,
			std::size_t const sample_ploidy_sum
		):
			task(),
			sequence_writer_task_delegate(),
			m_subgraph_semaphore(dispatch_semaphore_create(2 * hw_concurrency)), // FIXME: some other value?
			m_write_semaphore(dispatch_semaphore_create(2 * hw_concurrency)), // FIXME: some other value?
			m_path_permutation(generate_config.generated_path_count),
			m_reader(std::move(reader)),
			m_delegate(&delegate),
			m_generate_config(&generate_config),
			m_logger(&logger),
			m_preprocessing_result(&result),
			m_vcf_input_handle(&vcf_input_handle),
			m_record_count(record_count),
			m_sample_ploidy_sum(sample_ploidy_sum)
		{
			m_logger->status_logger.set_delegate(m_progress_counter);
			std::iota(m_path_permutation.begin(), m_path_permutation.end(), 0);
		}
		
		// task
		virtual void execute() override;
		
		// merge_subgraph_paths_task_delegate
		virtual void task_did_finish(
			merge_subgraph_paths_task &task,
			std::vector <reduced_subgraph::path_index> &&matchings,
			weight_type const matching_weight
		) override;
		virtual void task_did_calculate_edge_weight(merge_subgraph_paths_task &task) override { m_progress_counter.merge_subgraph_paths_task_did_calculate_edge_weight(); }
		
		// read_subgraph_variants_task_delegate
		virtual void task_did_finish(read_subgraph_variants_task &task, reduced_subgraph &&rsg) override;	
		virtual void task_did_handle_variant(read_subgraph_variants_task &task, variant const &variant) override { m_progress_counter.read_subgraph_variants_task_did_handle_variant(); }
		
		// sequence_writer_task_delegate
		virtual void handled_all_haplotypes(sequence_writer_task &task) override;
		virtual void task_did_finish(sequence_writer_task &task) override;
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) override;
			
	protected:
		std::size_t count_variants_in_range(
			std::size_t current_lineno,
			std::size_t const var_lineno,
			std::vector <std::size_t> &skipped_lines
		);
		
		void add_graph_range(
			std::size_t const start_lineno,
			std::size_t const end_lineno,
			std::size_t	const range_start_offset,
			std::size_t	const range_length,
			std::size_t const variant_count,
			std::vector <std::size_t> const &skipped_lines,
			std::vector <graph_range> &graph_ranges
		);
		
		void init_read_subgraph_variants_task(
			read_subgraph_variants_task &task,
			std::size_t const task_idx,
			graph_range &&range
		);
			
		void start_merge_task(std::size_t const lhs_idx, std::size_t const rhs_idx);
		auto create_haplotype(
			std::size_t const sample_no,
			std::size_t const ploidy
		) -> haplotype_map_type::iterator;
		void prepare_haplotypes();
	};
}


#endif
