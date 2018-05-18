/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GENERATE_SEQUENCES_CONTEXT_HH
#define VCF2MULTIALIGN_GENERATE_SEQUENCES_CONTEXT_HH

#include <vcf2multialign/generate_context_base.hh>
#include <vcf2multialign/variant_handler.hh>


namespace vcf2multialign {
	
	struct haplotype_os : public haplotype
	{
		libbio::file_ostream output_stream;
		
		void fill_with_dashes(std::size_t fill_amt)
		{
			std::ostream_iterator <char> it(output_stream);
			std::fill_n(it, fill_amt, '-');
		}
		
		void write(char const *output_start, std::size_t const length)
		{
			output_stream.write(output_start, length);
		}
		
		void append(std::string const &alt)
		{
			output_stream << alt;
		}
	};
	
	
	class generate_sequences_context final : public generate_context_base
	{
	protected:
		typedef variant_handler <haplotype_os, generate_sequences_context>	variant_handler_type;
		typedef typename variant_handler_type::haplotype_map_type			haplotype_map;
		
		friend variant_handler_type;
		
	protected:
		variant_handler_type												m_variant_handler;
		haplotype_map														m_haplotypes;
		libbio::vcf_reader::sample_name_map::const_iterator					m_sample_names_it{};
		libbio::vcf_reader::sample_name_map::const_iterator					m_sample_names_end{};
		std::size_t															m_chunk_size{0};
		std::size_t															m_current_round{0};
		std::size_t															m_total_rounds{0};
		bool																m_should_overwrite_files{false};
		
	public:
		generate_sequences_context(
			libbio::dispatch_ptr <dispatch_queue_t> &&main_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &&parsing_queue,
			char const *null_allele_seq,
			std::size_t const chunk_size,
			bool const should_overwrite_files
		):
			generate_context_base(
				std::move(main_queue),
				std::move(parsing_queue),
				null_allele_seq,
				should_overwrite_files
			),
			m_should_overwrite_files(should_overwrite_files)
		{
		}
		
		void load_and_generate(
			char const *reference_fname,
			char const *variants_fname,
			char const *out_reference_fname,
			char const *report_fname,
			sv_handling const sv_handling_method,
			bool const should_check_ref
		);
		
		
	protected:
		void update_haplotypes(char const *out_reference_fname);
		void generate_sequences(char const *out_reference_fname = nullptr);
		void finish_round();
		void variant_handler_did_finish(variant_handler_type &handler);
		void variant_handler_did_empty_overlap_stack(variant_handler_type &) {}
	};
}

#endif
