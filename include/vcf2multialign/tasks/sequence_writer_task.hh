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
		
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) = 0;
	};
	
	
	class sequence_writer_task :
		public parsing_task,
		public variant_stats,
		public sequence_writer_delegate <variant>
	{
	protected:
		sequence_writer <channel_ostream, variant>	m_sequence_writer;
		sequence_writer_task_delegate				*m_delegate{nullptr};
		alt_checker const							*m_alt_checker{nullptr};
		
	public:
		sequence_writer_task(
			sequence_writer_task_delegate &delegate,
			class status_logger &status_logger,
			class error_logger &error_logger,
			class vcf_reader const &vcf_reader,
			class alt_checker const &alt_checker,
			vector_type const &reference,
			std::string const &null_allele_seq
		):
			parsing_task(status_logger, error_logger, vcf_reader),
			m_sequence_writer(reference, null_allele_seq),
			m_delegate(&delegate),
			m_alt_checker(&alt_checker)
		{
			m_sequence_writer.set_delegate(*this);
		}
		
		void prepare(haplotype_map <channel_ostream> &haplotypes) { m_sequence_writer.prepare(haplotypes); }
		virtual void execute() override;
		
		// variant_stats
		virtual class error_logger &error_logger() override { return *m_error_logger; }
		virtual class status_logger &status_logger() override { return *m_status_logger; }
		
		// sequence_writer_delegate
		virtual std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const override { return m_alt_checker->valid_alts(lineno); }
		virtual bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const override { return m_alt_checker->is_valid_alt(lineno, alt_idx); }
		virtual void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		) override { m_delegate->enumerate_sample_genotypes(var, cb); }
		// Rest comes from variant_stats.
		
	};
}

#endif
