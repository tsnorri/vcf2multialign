/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_SEQUENCES_GENERATE_GRAPH_CONTEXT_HH
#define VCF2MULTIALIGN_SEQUENCES_GENERATE_GRAPH_CONTEXT_HH

#include <vcf2multialign/generate_context_base.hh>
#include <vcf2multialign/sequences/graph_writer.hh>
#include <vcf2multialign/variant_handler.hh>


namespace vcf2multialign {
	
	struct haplotype_vs : public haplotype
	{
		std::vector <uint8_t> output_sequence;
		
		void fill_with_dashes(std::size_t fill_amt)
		{
			output_sequence.resize(output_sequence.size() + fill_amt, '-');
		}
		
		void write(char const *output_start, std::size_t const length)
		{
			auto const pos(output_sequence.size());
			output_sequence.resize(pos + length, '\0');
			std::copy_n(output_start, length, output_sequence.begin() + pos);
		}
		
		void append(std::string const &alt)
		{
			auto const pos(output_sequence.size());
			auto const length(alt.size());
			output_sequence.resize(pos + length, '\0');
			std::copy_n(alt.data(), length, output_sequence.begin() + pos);
		}
	};
	
	std::ostream &operator<<(std::ostream &os, haplotype_vs const &haplotype);
	
	
	class generate_graph_context final : public generate_context_base
	{
	protected:
		typedef variant_handler <haplotype_vs, generate_graph_context>	variant_handler_type;
		typedef typename variant_handler_type::haplotype_map_type		haplotype_map;
		typedef graph_writer <generate_graph_context>					graph_writer_type;
		
		friend variant_handler_type;
		friend graph_writer_type;
		
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>		m_output_queue{};
		libbio::dispatch_ptr <dispatch_semaphore_t>	m_output_sema{};
		libbio::file_ostream						m_output_stream;
		variant_handler_type						m_variant_handler;
		haplotype_map 								m_haplotypes;	// FIXME: check if m_haplotypes could be a vector instead.
		std::vector <std::vector <uint8_t>>			m_buffers;
		graph_writer_type							m_graph_writer;
		
	public:
		generate_graph_context(
			libbio::dispatch_ptr <dispatch_queue_t> &&main_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &&parsing_queue,
			libbio::dispatch_ptr <dispatch_queue_t> &&output_queue,
			char const *chromosome_name,
			char const *null_allele_seq,
			bool const should_overwrite_files
		):
			generate_context_base(
				std::move(main_queue),
				std::move(parsing_queue),
				chromosome_name,
				null_allele_seq,
				should_overwrite_files
			),
			m_output_queue(std::move(output_queue)),
			m_output_sema(dispatch_semaphore_create(1), false)
		{
		}
		
	public:
		void load_and_generate(
			char const *reference_fname,
			char const *ref_seq_name,
			char const *variants_fname,
			char const *report_fname,
			char const *output_fname,
			sv_handling const sv_handling_method,
			bool const should_check_ref
		);
		
	protected:
		void open_files(
			char const *reference_fname,
			char const *ref_seq_name,
			char const *variants_fname,
			char const *report_fname,
			char const *output_fname
		);
		
		void prepare_haplotypes();
		void swap_buffers_and_generate_graph();
		void graph_writer_did_process_segment(graph_writer_type &);
		void variant_handler_did_empty_overlap_stack(variant_handler_type &);
		void variant_handler_did_finish(variant_handler_type &handler);
	};
}

#endif
