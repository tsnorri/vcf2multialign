/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_READ_SUBGRAPH_VARIANTS_TASK_HH
#define VCF2MULTIALIGN_TASKS_READ_SUBGRAPH_VARIANTS_TASK_HH

#include <sdsl/bits.hpp>
#include <vcf2multialign/generate_configuration.hh>
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
		typedef reduced_subgraph::sequence_type		sequence_type;
		typedef std::map <sample_id, sequence_type>	sequence_map;
		
	protected:
		read_subgraph_variants_task_delegate		*m_delegate{nullptr};
		sequence_map								m_sequences_by_sample;		// Used only in the first phase (before finish()).
		
		// End co-ordinate of the REF sequence in each sample.
		// Used only in the first phase (before finish()).
		std::vector <std::size_t>					m_endpoints;
		graph_range									m_subgraph_range;
		vcf_mmap_input								m_vcf_input;				// FIXME: move to parsing_task_vh?
		std::size_t									m_task_idx{0};
		std::size_t									m_generated_path_count{0};
		std::size_t									m_start_lineno{0};
		std::size_t									m_variant_count{0};
		std::size_t									m_alt_field_width{0};
		
	public:
		read_subgraph_variants_task() = default;
		
		read_subgraph_variants_task(
			read_subgraph_variants_task_delegate &delegate,
			generate_configuration const &generate_config,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,	// Needs to be serial.
			struct logger &logger,
			class vcf_reader const &vcf_reader,
			mmap_handle const &vcf_input_handle,
			alt_checker const &checker,
			vector_type const &reference,
			variant_set const &skipped_variants,
			graph_range &&range,
			std::size_t const task_idx, 
			std::size_t const sample_ploidy_sum
		):
			parsing_task_vh(
				worker_queue,
				logger,
				vcf_reader,
				checker,
				reference,
				generate_config.sv_handling_method,
				skipped_variants
			),
			m_delegate(&delegate),
			m_endpoints(sample_ploidy_sum, 0),
			m_subgraph_range(std::move(range)),
			m_vcf_input(vcf_input_handle),
			m_task_idx(task_idx),
			m_generated_path_count(generate_config.generated_path_count),
			m_alt_field_width(1 << sdsl::bits::hi(checker.max_alt_field_size()))
		{
			parsing_task_vh::vcf_reader().set_input(m_vcf_input);
		}
		
		std::size_t task_idx() const { return m_task_idx; }
		virtual void execute() override;
		
		// variant_handler_delegate
		virtual void prepare(class vcf_reader &reader) override;
		void handle_variant(variant &var) override;
		void finish() override;
		
		// variant_handler_container
		virtual void finish_copy_or_move() override
		{
			parsing_task_vh::finish_copy_or_move();
			vcf_reader().set_input(m_vcf_input);
		}
		
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
