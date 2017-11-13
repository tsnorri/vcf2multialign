/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_SEQUENCE_WRITER_TASK_HH
#define VCF2MULTIALIGN_TASKS_SEQUENCE_WRITER_TASK_HH

#include <vcf2multialign/sequence_writer.hh>
#include <vcf2multialign/tasks/parsing_task.hh>
#include <vcf2multialign/variant_stats.hh>


namespace vcf2multialign {
	
	class sequence_writer_task;
	
	
	struct sequence_writer_task_delegate
	{
		virtual ~sequence_writer_task_delegate() {}
		
		virtual void task_did_finish(sequence_writer_task &task) = 0;
		virtual void handled_all_haplotypes(sequence_writer_task &task) = 0;
		
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) = 0;
	};
	
	
	class sequence_writer_task :
		public parsing_task_vh,
		public variant_stats,
		public sequence_writer_delegate <variant>
	{
	protected:
		typedef sequence_writer <channel_ostream, variant> sequence_writer_type;
			
	protected:
		sequence_writer_type						m_sequence_writer;
		sequence_writer_task_delegate				*m_delegate{nullptr};
		
	public:
		sequence_writer_task(
			sequence_writer_task_delegate &delegate,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,
			class status_logger &status_logger,
			class error_logger &error_logger,
			class vcf_reader const &vcf_reader,
			class alt_checker const &alt_checker,
			variant_set const &skipped_variants,
			vector_type const &reference,
			std::string const &null_allele_seq,
			sv_handling const sv_handling_method
		):
			parsing_task_vh(
				worker_queue,
				status_logger,
				error_logger,
				vcf_reader,
				alt_checker,
				reference,
				sv_handling_method,
				skipped_variants
			),
			m_sequence_writer(reference, null_allele_seq),
			m_delegate(&delegate)
		{
			m_sequence_writer.set_delegate(*this);
		}
		
		void prepare(haplotype_map <channel_ostream> &haplotypes) { m_sequence_writer.prepare(haplotypes); }
		virtual void execute() override;
		
		// variant_stats
		virtual class error_logger &error_logger() override { return *m_error_logger; }
		virtual class status_logger &status_logger() override { return *m_status_logger; }
		
		// variant_handler_delegate
		virtual void prepare(class vcf_reader &reader) override;
		void handle_variant(variant &var) override;
		void finish() override;
		
		// sequence_writer_delegate
		virtual std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const override { return m_alt_checker->valid_alts(lineno); }
		virtual bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const override { return m_alt_checker->is_valid_alt(lineno, alt_idx); }
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) override { m_delegate->enumerate_sample_genotypes(var, cb); }
		virtual void handled_all_haplotypes() override { m_delegate->handled_all_haplotypes(*this); }
		// Rest comes from variant_stats.
	};
}

#endif
