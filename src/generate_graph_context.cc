/*
 * Copyright (c) 2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/generate_graph_context.hh>


namespace lb = libbio;


namespace vcf2multialign {
	
	void generate_graph_context::open_files(
		char const *reference_fname,
		char const *variants_fname,
		char const *report_fname,
		char const *output_fname
	)
	{
		generate_context_base::open_files(reference_fname, variants_fname, report_fname);
		
		lb::open_file_for_writing(output_fname, m_output_stream, m_should_overwrite_files);
	}
	
	
	void generate_graph_context::swap_buffers_and_generate_graph()
	{
		// Make sure that m_graph_writer has finished with the previous set of buffers.
		dispatch_semaphore_wait(*m_output_sema, DISPATCH_TIME_FOREVER);

		// Swap buffers between haplotypes and m_buffers.
		// Also reset the range pointers in the new haplotype buffers.
		std::size_t i(0);
		for (auto &kv : m_haplotypes)
		{
			auto &haplotype_vec(kv.second);
			for (auto &haplotype : haplotype_vec)
			{
				assert(i < m_buffers.size());
				auto &buffer(m_buffers[i]);
				
				using std::swap;
				swap(haplotype.output_sequence, buffer);
				
				haplotype.output_sequence.clear();
				++i;
			}
		}
		
		// Process in background thread.
		lb::dispatch_async_fn(*m_output_queue, [this](){
			m_graph_writer.process_segment(m_buffers);
		});
	}
	
	
	void generate_graph_context::graph_writer_did_process_segment(graph_writer_type &)
	{
		// Done with the previous segment, buffers are writable again.
		dispatch_semaphore_signal(*m_output_sema);
	}
	
	
	void generate_graph_context::variant_handler_did_process_overlap_stack(variant_handler_type &handler)
	{
		if (1 == handler.m_overlap_stack.size())
			swap_buffers_and_generate_graph();
	}
	
	
	void generate_graph_context::variant_handler_did_finish(variant_handler_type &handler)
	{
		// Handle the last segment.
		swap_buffers_and_generate_graph();
		
		// Finalize the graph.
		lb::dispatch_caller caller(&m_graph_writer);
		caller.template sync <&graph_writer_type::finish>(*m_output_queue);
		
		// Everything done.
		finish();
	}
	
	
	void generate_graph_context::prepare_haplotypes()
	{
		// Initialize the contents of m_haplotypes.
		std::size_t count(0);
		auto const &sample_names(m_vcf_reader.sample_names());
		for (auto const &kv : sample_names)
		{
			auto const &sample_name(kv.first);
			auto const sample_no(kv.second);
			auto const current_ploidy(m_ploidy.find(sample_no)->second);
			
			auto const did_emplace(m_haplotypes.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(sample_no),
				std::forward_as_tuple(current_ploidy)
			).second);
			lb::always_assert(did_emplace);
			count += current_ploidy;
		}
		
		{
			lb::always_assert(
				m_haplotypes.cend() == m_haplotypes.find(REF_SAMPLE_NUMBER),
				"REF_SAMPLE_NUMBER already in use"
			);
			
			auto const did_emplace(m_haplotypes.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(REF_SAMPLE_NUMBER),
				std::forward_as_tuple(1)
			).second);
			lb::always_assert(did_emplace);
			++count;
		}
		
		m_buffers.resize(count);
		m_graph_writer.init(count);
	}
	
	
	void generate_graph_context::load_and_generate(
		char const *reference_fname,
		char const *variants_fname,
		char const *report_fname,
		char const *output_fname,
		sv_handling const sv_handling_method,
		bool const should_check_ref
	)
	{
		// Open the files.
		std::cerr << "Opening files…" << std::endl;
		open_files(reference_fname, variants_fname, report_fname, output_fname);
		
		generate_context_base::load_and_generate(
			sv_handling_method,
			should_check_ref
		);
		
		// Replace the placeholder variant_handler.
		{
			variant_handler_type temp(
				m_main_queue,
				m_parsing_queue,
				m_vcf_reader,
				m_reference,
				sv_handling_method,
				m_skipped_variants,
				m_null_allele_seq,
				m_error_logger,
				*this
			);
		
			m_variant_handler = std::move(temp);
			m_variant_handler.get_variant_buffer().set_delegate(m_variant_handler);
		}
		
		// Replace the placeholder graph_writer.
		{
			graph_writer_type temp(*this, m_output_stream);
			m_graph_writer = std::move(temp);
		}
		
		// Create the required buffers.
		prepare_haplotypes();
		
		std::cerr << "Generating variant graph…" << std::endl;
		m_start_time = std::chrono::system_clock::now();
		m_variant_handler.process_variants(m_haplotypes);
	}
}
