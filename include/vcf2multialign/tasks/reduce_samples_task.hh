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
		public sequence_writer_task_delegate,
		public status_logger_delegate
	{
	public:
		typedef haplotype_map <channel_ostream>		haplotype_map_type;
		typedef std::vector <reduced_subgraph>		subgraph_vector;
		typedef std::vector <
			std::vector <
				reduced_subgraph::path_index
			>
		>											path_matchings_type;
		
	protected:
		dispatch_ptr <dispatch_semaphore_t>			m_semaphore{};	// FIXME: make use of this.
		dispatch_ptr <dispatch_semaphore_t>			m_write_semaphore{};

		std::mutex									m_subgraph_mutex;
		subgraph_vector								m_subgraphs;
		std::vector <bool>							m_subgraph_bitmap;
		path_matchings_type							m_path_matchings;
		std::vector <std::size_t>					m_path_permutation;
		haplotype_map_type							m_haplotypes;
		subgraph_vector::const_iterator				m_subgraph_iterator{};
		vcf_reader									m_reader;
		
		reduce_samples_task_delegate				*m_delegate{nullptr};
		status_logger								*m_status_logger{nullptr};
		error_logger								*m_error_logger{nullptr};
		vector_type const							*m_reference{nullptr};
		std::string const							*m_null_allele_seq{nullptr};
		alt_checker const							*m_alt_checker{nullptr};
		subgraph_map const							*m_subgraph_starting_points{nullptr};
		variant_set const							*m_skipped_variants{nullptr};
		boost::optional <std::string> const			*m_out_reference_fname{nullptr};
		sv_handling									m_sv_handling_method{};
		std::atomic_size_t							m_remaining_merge_tasks{0};
		std::size_t									m_record_count{0};
		std::size_t									m_generated_path_count{0};
		std::size_t									m_sample_ploidy_sum{0};
		std::size_t									m_min_path_length{0};
		bool										m_should_overwrite_files{false};
		
	public:
		reduce_samples_task() = default;
		
		reduce_samples_task(
			reduce_samples_task_delegate &delegate,
			status_logger &status_logger,
			error_logger &error_logger,
			std::size_t const hw_concurrency,
			vcf_reader &&reader,
			vector_type const &reference,
			std::string const &null_allele_seq,
			alt_checker const &checker,
			subgraph_map const &subgraph_starting_points,
			variant_set const &skipped_variants,
			boost::optional <std::string> const &out_reference_fname,
			sv_handling const sv_handling_method,
			std::size_t const record_count,
			std::size_t const generated_path_count,
			std::size_t const sample_ploidy_sum,
			std::size_t const min_path_length,
			bool const should_overwrite_files
		):
			task(),
			sequence_writer_task_delegate(),
			status_logger_delegate(),
			m_semaphore(dispatch_semaphore_create(2 * hw_concurrency)), // FIXME: some other value?
			m_write_semaphore(dispatch_semaphore_create(2 * hw_concurrency)), // FIXME: arbitrary value.
			m_path_permutation(generated_path_count),
			m_reader(std::move(reader)),
			m_delegate(&delegate),
			m_status_logger(&status_logger),
			m_error_logger(&error_logger),
			m_reference(&reference),
			m_null_allele_seq(&null_allele_seq),
			m_alt_checker(&checker),
			m_subgraph_starting_points(&subgraph_starting_points),
			m_skipped_variants(&skipped_variants),
			m_out_reference_fname(&out_reference_fname),
			m_sv_handling_method(sv_handling_method),
			m_record_count(record_count),
			m_generated_path_count(generated_path_count),
			m_sample_ploidy_sum(sample_ploidy_sum),
			m_min_path_length(min_path_length),
			m_should_overwrite_files(should_overwrite_files)
		{
			m_status_logger->set_delegate(*this);
			std::iota(m_path_permutation.begin(), m_path_permutation.end(), 0);
		}
		
		// task
		virtual void execute() override;
		
		// status_logger_delegate
		virtual std::size_t record_count() const override { return m_record_count; }
		virtual std::size_t current_record() const override { return 0; } // FIXME return the correct value.
		
		// merge_subgraph_paths_task_delegate
		virtual void task_did_finish(
			merge_subgraph_paths_task &task,
			std::vector <reduced_subgraph::path_index> &&matchings
		) override;
		
		// read_subgraph_variants_task_delegate
		virtual void task_did_finish(read_subgraph_variants_task &task, reduced_subgraph &&rsg) override;	
		
		// sequence_writer_task_delegate
		virtual void task_did_finish(sequence_writer_task &task) override;
		
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) override;
			
	protected:
		void init_read_subgraph_variants_task(
			read_subgraph_variants_task &task,
			std::size_t const task_idx,
			char const *buffer_start,
			char const *buffer_end,
			std::size_t const start_lineno,
			std::size_t const variant_count
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
