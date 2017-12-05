/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_ALL_HAPLOTYPES_TASK_HH
#define VCF2MULTIALIGN_TASKS_ALL_HAPLOTYPES_TASK_HH

#include <vcf2multialign/generate_configuration.hh>
#include <vcf2multialign/logger.hh>
#include <vcf2multialign/preprocessing_result.hh>
#include <vcf2multialign/tasks/sequence_writer_task.hh>


namespace vcf2multialign {
	
	class all_haplotypes_task;
	
	
	struct all_haplotypes_task_delegate
	{
		virtual ~all_haplotypes_task_delegate() {}
		virtual void task_did_finish(all_haplotypes_task &task) = 0;
		virtual void store_and_execute(std::unique_ptr <task> &&task) = 0;
		virtual void remove_task(task &task) = 0;
	};
	
	
	class all_haplotypes_task :
		public task,
		public sequence_writer_task_delegate,
		public status_logger_delegate
	{
	protected:
		typedef haplotype_map <channel_ostream>				haplotype_map_type;
		
	protected:
		dispatch_ptr <dispatch_semaphore_t>					m_write_semaphore;

		all_haplotypes_task_delegate						*m_delegate{};
		logger												*m_logger{};
		generate_configuration const						*m_generate_config{};
		preprocessing_result const							*m_preprocessing_result{};
		
		vcf_reader											m_vcf_reader;
		haplotype_map_type									m_haplotypes;
		
		vcf_reader::sample_name_map::const_iterator			m_sample_names_it{};
		vcf_reader::sample_name_map::const_iterator			m_sample_names_end{};
		
		// FIXME: change all instances of std::chrono::system_clock to steady_clock when measuring elapsed time.
		std::chrono::time_point <std::chrono::system_clock>	m_start_time{};
		std::chrono::time_point <std::chrono::system_clock>	m_round_start_time{};
		
		std::size_t											m_record_count{};
		std::size_t											m_current_round{};
		std::size_t											m_total_rounds{};
		
	public:
		all_haplotypes_task(
			all_haplotypes_task_delegate &delegate,
			generate_configuration const &generate_config,
			struct logger &logger,
			class vcf_reader &&vcf_reader,
			preprocessing_result const &result,
			std::size_t const hw_concurrency,
			std::size_t const record_count
		):
			m_write_semaphore(dispatch_semaphore_create(2 * hw_concurrency)), // FIXME: some other value?
			m_delegate(&delegate),
			m_logger(&logger),
			m_generate_config(&generate_config),
			m_preprocessing_result(&result),
			m_vcf_reader(std::move(vcf_reader)),
			m_record_count(record_count)
		{
			m_logger->status_logger.set_delegate(*this);
		}
		
		// task
		virtual void execute() override;
		
		// status_logger_delegate
		virtual std::size_t step_count() const override { return m_record_count; }
		virtual std::size_t current_step() const override { return vcf_reader().counter_value(); }
		
		// sequence_writer_task_delegate
		virtual void task_did_finish(sequence_writer_task &task) override;
		virtual void handled_all_haplotypes(sequence_writer_task &task) override {}
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) override;
		
	protected:
		haplotype_map_type::iterator create_haplotype(
			std::size_t const sample_no,
			std::size_t const ploidy
		);
	
		haplotype_map_type::iterator find_or_create_haplotype(
			std::size_t const sample_no,
			std::size_t const ploidy
		);
			
		void update_haplotypes_with_ref();
		void update_haplotypes(bool const output_reference);
		
		bool prepare_next_round(sequence_writer_task &task, bool const output_reference);
		void do_finish(sequence_writer_task &task);
	};
}

#endif
