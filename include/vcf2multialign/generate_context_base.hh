/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GENERATE_CONTEXT_BASE_HH
#define VCF2MULTIALIGN_GENERATE_CONTEXT_BASE_HH

#include <chrono>
#include <map>
#include <libbio/dispatch.hh>
#include <libbio/mmap_handle.hh>
#include <libbio/vcf/vcf_reader.hh>
#include <string>
#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class generate_context_base
	{
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>				m_main_queue{};
		libbio::dispatch_ptr <dispatch_queue_t>				m_parsing_queue{};
		vector_type											m_reference;
		libbio::mmap_handle <char>							m_vcf_handle;
		libbio::vcf_mmap_input								m_vcf_input;
		libbio::vcf_reader									m_vcf_reader;
		
		std::chrono::time_point <std::chrono::system_clock>	m_start_time{};
		std::chrono::time_point <std::chrono::system_clock>	m_round_start_time{};
		
		error_logger										m_error_logger;
		
		ploidy_map											m_ploidy;
		variant_set											m_skipped_variants;

		std::string											m_chromosome_name;
		std::string											m_null_allele_seq;
		bool												m_should_overwrite_files{false};
		
	public:
		generate_context_base(
			libbio::dispatch_ptr <dispatch_queue_t> &&main_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &&parsing_queue,
			char const *chromosome_name,
			char const *null_allele_seq,
			bool const should_overwrite_files
		):
			m_main_queue(std::move(main_queue)),
			m_parsing_queue(std::move(parsing_queue)),
			m_vcf_input(m_vcf_handle),	// XXX makes copying and moving this class non-trivial.
			m_chromosome_name(chromosome_name ?: ""),
			m_null_allele_seq(null_allele_seq),
			m_should_overwrite_files(should_overwrite_files)
		{
		}
		
		generate_context_base(generate_context_base const &) = delete;
		generate_context_base(generate_context_base &&) = delete;
		
		virtual ~generate_context_base() {}
		
	protected:
		void finish();
		void cleanup() { delete this; }
		void open_files(
			char const *reference_fname,
			char const *ref_seq_name,
			char const *variants_fname,
			char const *report_fname
		);
		void check_ploidy();
		void check_ref();
		void load_and_generate(bool const should_check_ref);
	};
}

#endif
