/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_ALL_HAPLOTYPES_TASK_HH
#define VCF2MULTIALIGN_TASKS_ALL_HAPLOTYPES_TASK_HH

#include <vcf2multialign/sequence_writer.hh>
#include <vcf2multialign/tasks/parsing_task.hh>
#include <vcf2multialign/variant_stats.hh>


namespace vcf2multialign {
	
	class all_haplotypes_task;
	
	
	struct all_haplotypes_task_delegate
	{
		virtual void task_did_finish(all_haplotypes_task &task) = 0;
	};
	
	
	class all_haplotypes_task :
		public parsing_task_vh,
		public variant_stats,
		public status_logger_delegate,
		public virtual sequence_writer_delegate
	{
	protected:
		all_haplotypes_task_delegate						*m_delegate{nullptr};
		
		ploidy_map const									*m_ploidy{nullptr};
		variant_set const									*m_skipped_variants{nullptr};
		boost::optional <std::string> const					*m_out_reference_fname{nullptr};
		
		haplotype_map										m_haplotypes;
		sequence_writer										m_sequence_writer;
		
		vcf_reader::sample_name_map::const_iterator			m_sample_names_it{};
		vcf_reader::sample_name_map::const_iterator			m_sample_names_end{};
		
		std::chrono::time_point <std::chrono::system_clock>	m_start_time{};
		std::chrono::time_point <std::chrono::system_clock>	m_round_start_time{};
		
		std::size_t											m_record_count{0};
		std::size_t											m_current_round{0};
		std::size_t											m_total_rounds{0};
		std::size_t											m_chunk_size{0};
		bool												m_should_overwrite_files{false};
		
	public:
		all_haplotypes_task() = delete;
		all_haplotypes_task(all_haplotypes_task const &) = delete;
		all_haplotypes_task(all_haplotypes_task &&) = delete;
		all_haplotypes_task &operator=(all_haplotypes_task const &) = delete;
		all_haplotypes_task &operator=(all_haplotypes_task &&) = delete;
		
		all_haplotypes_task(
			all_haplotypes_task_delegate &delegate,
			class status_logger &status_logger,
			class error_logger &error_logger,
			class variant_handler &&variant_handler,
			class vcf_reader &&vcf_reader,
			vector_type const &reference,
			std::string const &null_allele_seq,
			ploidy_map const &ploidy,
			variant_set const &skipped_variants,
			boost::optional <std::string> const &out_reference_fname,
			std::size_t record_count,
			std::size_t chunk_size,
			bool should_overwrite_files
		):
			parsing_task_vh(status_logger, error_logger, std::move(variant_handler), std::move(vcf_reader)),
			m_delegate(&delegate),
			m_ploidy(&ploidy),
			m_skipped_variants(&skipped_variants),
			m_out_reference_fname(&out_reference_fname),
			m_sequence_writer(reference, null_allele_seq),
			m_chunk_size(chunk_size),
			m_should_overwrite_files(should_overwrite_files)
		{
			m_status_logger->set_delegate(*this);
			m_sequence_writer.set_delegate(*this);
		}
		
		// variant_handler_delegate (from parsing_task_vh_base)
		virtual void handle_variant(variant &var) override;
		virtual void finish() override;
		
		// parsing_task
		virtual void execute() override;
		
		// variant_stats
		virtual class error_logger &error_logger() override { return *m_error_logger; }
		virtual class status_logger &status_logger() override { return *m_status_logger; }
		
		// status_logger_delegate
		virtual std::size_t record_count() const override { return m_record_count; }
		virtual std::size_t current_record() const override { return m_vcf_reader.counter_value(); }
		
		// sequence_writer_delegate
		virtual std::set <std::size_t> const &valid_alts(variant &var) const override { return m_variant_handler.valid_alts(); }
		virtual bool is_valid_alt(std::size_t const alt_idx) const override { return m_variant_handler.is_valid_alt(alt_idx); }
		virtual void enumerate_genotype(
			variant &var,
			std::size_t const sample_no,
			std::function <void(uint8_t, std::size_t, bool)> const &cb
		) override { m_variant_handler.enumerate_genotype(var, sample_no, cb); }
		// Rest comes from variant_stats.
		
	protected:
		haplotype_map::iterator create_haplotype(
			std::size_t const sample_no,
			std::size_t const ploidy
		);
	
		haplotype_map::iterator find_or_create_haplotype(
			std::size_t const sample_no,
			std::size_t const ploidy
		);
			
		void update_haplotypes_with_ref();
		void update_haplotypes(bool const output_reference);
		void generate_sequences(bool const output_reference = false);
	};
}

#endif
