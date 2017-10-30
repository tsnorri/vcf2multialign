/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_READ_SUBGRAPH_VARIANTS_TASK_HH
#define VCF2MULTIALIGN_TASKS_READ_SUBGRAPH_VARIANTS_TASK_HH

#include <vcf2multialign/reduced_subgraph.hh>
#include <vcf2multialign/tasks/parsing_task.hh>


namespace vcf2multialign {
	
	class read_subgraph_variants_task;
	
	
	struct read_subgraph_variants_task_delegate
	{
		virtual ~read_subgraph_variants_task_delegate() {}
		virtual void task_did_finish(read_subgraph_variants_task &task, reduced_subgraph &&rsg) = 0;
		virtual void task_did_handle_variant(read_subgraph_variants_task &task, variant const &variant) = 0;
	};
	
	
	class read_subgraph_variants_task :
		public parsing_task_vh
	{
	protected:
		read_subgraph_variants_task_delegate	*m_delegate{nullptr};
		reduced_subgraph::sequence_map			m_sequences_by_sample;		// Used only in the first phase (before finish()).
		
		// End co-ordinate of the REF sequence in each sample.
		// Used only in the first phase (before finish()).
		std::vector <std::size_t>				m_endpoints;
		std::size_t								m_task_idx{0};
		std::size_t								m_generated_path_count{0};
		std::size_t								m_start_lineno{0};
		std::size_t								m_variant_count{0};
		
	public:
		read_subgraph_variants_task() = default;
		
		read_subgraph_variants_task(
			read_subgraph_variants_task_delegate &delegate,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,	// Needs to be serial.
			status_logger &status_logger,
			error_logger &error_logger,
			class vcf_reader const &vcf_reader,
			alt_checker const &checker,
			vector_type const &reference,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants,
			boost::optional <std::string> const &out_reference_fname,
			std::size_t const task_idx, 
			std::size_t const generated_path_count,
			std::size_t const sample_ploidy_sum,
			std::size_t const start_lineno,
			std::size_t const variant_count
		):
			parsing_task_vh(
				worker_queue,
				status_logger,
				error_logger,
				vcf_reader,
				checker,
				reference,
				sv_handling_method,
				skipped_variants
			),
			m_delegate(&delegate),
			m_endpoints(sample_ploidy_sum, 0),
			m_task_idx(task_idx),
			m_generated_path_count(generated_path_count),
			m_start_lineno(start_lineno),
			m_variant_count(variant_count)
		{
		}
		
		std::size_t task_idx() const { return m_task_idx; }
		virtual void execute() override;
		
		// variant_handler_delegate
		void handle_variant(variant &var) override;
		void finish() override;
		
	protected:
		void invert_sequences_by_sample(
			reduced_subgraph::sequence_vec &sequences,
			reduced_subgraph::sample_bimap &samples_by_sequence_idx
		);
			
		void split_sequences_to_paths(
			std::size_t const original_seq_count,
			reduced_subgraph::sample_bimap const &samples_by_sequence_idx,
			reduced_subgraph::path_map &generated_paths,
			reduced_subgraph::path_eq_map &generated_paths_eq
		);
	};
}


#endif
