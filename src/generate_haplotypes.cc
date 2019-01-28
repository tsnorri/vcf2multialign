/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/generate_graph_context.hh>
#include <vcf2multialign/generate_haplotypes.hh>
#include <vcf2multialign/generate_sequences_context.hh>


namespace lb = libbio;


namespace vcf2multialign {
	
	void generate_haplotypes(
		char const *reference_fname,
		char const *variants_fname,
		char const *ref_seq_name,
		output_type const ot,
		char const *out_fname,
		char const *out_reference_fname,
		char const *report_fname,
		char const *null_allele_seq,
		std::size_t const chunk_size,
		sv_handling const sv_handling_method,
		bool const should_overwrite_files,
		bool const should_check_ref
	)
	{
		lb::dispatch_ptr <dispatch_queue_t> main_queue(dispatch_get_main_queue(), true);
		lb::dispatch_ptr <dispatch_queue_t> parsing_queue(
			dispatch_queue_create("fi.iki.tsnorri.vcf2multialign.parsing_queue", DISPATCH_QUEUE_SERIAL),
			false
		);
		
		// generate_context_base subclasses need to be allocated on the heap because later dispatch_main is called.
		// The class deallocates itself in cleanup().
		switch (ot)
		{
			case output_type::SEQUENCE_FILES:
			{
				auto *ctx(new generate_sequences_context(
					std::move(main_queue),
					std::move(parsing_queue),
					null_allele_seq,
					chunk_size,
					should_overwrite_files
				));
		
				ctx->load_and_generate(
					reference_fname,
					variants_fname,
					ref_seq_name,
					out_reference_fname,
					report_fname,
					sv_handling_method,
					should_check_ref
				);
				
				break;
			}
			
			case output_type::VARIANT_GRAPH:
			{
				lb::dispatch_ptr <dispatch_queue_t> output_queue(
					dispatch_queue_create("fi.iki.tsnorri.vcf2multialign.output_queue", DISPATCH_QUEUE_SERIAL),
					false
				);
				
				auto *ctx(new generate_graph_context(
					std::move(main_queue),
					std::move(parsing_queue),
					std::move(output_queue),
					null_allele_seq,
					should_overwrite_files
				));
				
				ctx->load_and_generate(
					reference_fname,
					variants_fname,
					ref_seq_name,
					report_fname,
					out_fname,
					sv_handling_method,
					should_check_ref
				);
				
				break;
			}
			
			default:
				libbio_fail("Unexpected output type.");
				break;
		}
	}
}
