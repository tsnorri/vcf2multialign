/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_ALL_HAPLOTYPES_TASK_HH
#define VCF2MULTIALIGN_TASKS_ALL_HAPLOTYPES_TASK_HH

#include <vcf2multialign/generate_configuration.hh>
#include <vcf2multialign/logger.hh>
#include <vcf2multialign/sequence_writer.hh>
#include <vcf2multialign/tasks/parsing_task.hh>
#include <vcf2multialign/variant_stats.hh>


namespace vcf2multialign {
	
	class all_haplotypes_task;
	
	
	struct all_haplotypes_task_delegate
	{
		virtual ~all_haplotypes_task_delegate() {}
		virtual void task_did_finish(all_haplotypes_task &task) = 0;
	};
	
	
	class all_haplotypes_task :
		public parsing_task_vh,
		public variant_stats,
		public status_logger_delegate,
		public sequence_writer_delegate <variant>
	{
	public:
		typedef haplotype_map <file_ostream>				haplotype_map_type;
		typedef sequence_writer <file_ostream, variant>		sequence_writer_type;
		
	protected:
		all_haplotypes_task_delegate						*m_delegate{nullptr};
		generate_configuration const						*m_generate_config{nullptr};
		
		ploidy_map const									*m_ploidy{nullptr};
		variant_set const									*m_skipped_variants{nullptr};
		
		haplotype_map_type									m_haplotypes;
		sequence_writer_type								m_sequence_writer;
		
		vcf_reader::sample_name_map::const_iterator			m_sample_names_it{};
		vcf_reader::sample_name_map::const_iterator			m_sample_names_end{};
		
		// FIXME: change all instances of std::chrono::system_clock to steady_clock when measuring elapsed time.
		std::chrono::time_point <std::chrono::system_clock>	m_start_time{};
		std::chrono::time_point <std::chrono::system_clock>	m_round_start_time{};
		
		std::size_t											m_record_count{0};
		std::size_t											m_current_round{0};
		std::size_t											m_total_rounds{0};
		
	public:
		all_haplotypes_task() = delete;
		all_haplotypes_task(all_haplotypes_task const &) = delete;
		all_haplotypes_task(all_haplotypes_task &&) = delete;
		all_haplotypes_task &operator=(all_haplotypes_task const &) = delete;
		all_haplotypes_task &operator=(all_haplotypes_task &&) = delete;
		
		all_haplotypes_task(
			all_haplotypes_task_delegate &delegate,
			generate_configuration const &generate_config,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,
			struct logger &logger,
			class vcf_reader &&vcf_reader,
			vector_type const &reference,
			alt_checker const &checker,
			ploidy_map const &ploidy,
			variant_set const &skipped_variants,
			std::size_t const record_count
		):
			parsing_task_vh(
				worker_queue,
				logger,
				std::move(vcf_reader),
				checker,
				reference,
				generate_config.sv_handling_method,
				skipped_variants
			),
			variant_stats(),
			m_delegate(&delegate),
			m_generate_config(&generate_config),
			m_ploidy(&ploidy),
			m_skipped_variants(&skipped_variants),
			m_sequence_writer(reference, generate_config.null_allele_seq),
			m_record_count(record_count)
		{
			m_logger->status_logger.set_delegate(*this);
			m_sequence_writer.set_delegate(*this);
		}
		
		// variant_handler_delegate (from parsing_task_vh_base)
		virtual void handle_variant(variant &var) override;
		virtual void finish() override;
		
		// task
		virtual void execute() override;
		
		// variant_stats
		virtual struct logger &logger() override { return *m_logger; }
		
		// status_logger_delegate
		virtual std::size_t step_count() const override { return m_record_count; }
		virtual std::size_t current_step() const override { return vcf_reader().counter_value(); }
		
		// sequence_writer_delegate
		virtual std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const override { return m_alt_checker->valid_alts(lineno); }
		virtual bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const override { return m_alt_checker->is_valid_alt(lineno, alt_idx); }
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) override { variant_handler().enumerate_sample_genotypes(var, cb); }
		// Rest comes from variant_stats.
		
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
		void generate_sequences(bool const output_reference = false);
	};
}

#endif
