/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_PREPARATION_TASK_HH
#define VCF2MULTIALIGN_TASKS_PREPARATION_TASK_HH

#include <vcf2multialign/alt_checker.hh>
#include <vcf2multialign/status_logger.hh>
#include <vcf2multialign/tasks/parsing_task.hh>


namespace vcf2multialign {
	
	class preparation_task;
	
	
	struct preparation_task_delegate
	{
		virtual void task_did_finish(preparation_task &task) = 0;
	};
	
	
	class preparation_task : public parsing_task, public status_logger_delegate
	{
	protected:
		preparation_task_delegate	*m_delegate{nullptr};
		vector_type const			*m_reference{nullptr};
		
		ploidy_map					m_ploidy;
		variant_set					m_skipped_variants;
		subgraph_map				m_subgraph_starting_points;
		alt_checker					m_alt_checker;
		
		std::size_t					m_record_count{0};
		sv_handling					m_sv_handling_method;
		bool						m_should_check_ref{false};
		
	public:
		preparation_task() = delete;
		
		preparation_task(
			preparation_task_delegate &delegate,
			status_logger &status_logger,
			error_logger &error_logger,
			vector_type const &reference,
			class vcf_reader &&vcf_reader,
			sv_handling const sv_handling_method,
			bool const should_check_ref
		):
			parsing_task(status_logger, error_logger, std::move(vcf_reader)),
			m_delegate(&delegate),
			m_reference(&reference),
			m_sv_handling_method(sv_handling_method),
			m_should_check_ref(should_check_ref)
		{
			m_status_logger->set_delegate(*this);
		}
		
		ploidy_map &ploidy_map() { return m_ploidy; }
		variant_set &skipped_variants() { return m_skipped_variants; }
		subgraph_map &subgraph_starting_points() { return m_subgraph_starting_points; }
		alt_checker &alt_checker() { return m_alt_checker; }
		
		// parsing_task
		virtual void execute() override;
		
		// status_logger_delegate
		virtual std::size_t record_count() const override { return m_record_count; }
		virtual std::size_t current_record() const override { return m_vcf_reader.counter_value(); }
		
	protected:
		void check_ploidy();
		void check_ref();
	};
}


#endif
