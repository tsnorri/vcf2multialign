/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/check_overlapping_non_nested_variants.hh>
#include <vcf2multialign/find_subgraph_starting_points.hh>
#include <vcf2multialign/tasks/preparation_task.hh>


namespace v2m = vcf2multialign;


namespace {
	bool compare_references(v2m::vector_type const &ref, std::string_view const &var_ref, std::size_t const var_pos, std::size_t /* out */ &idx);
}


namespace {
	bool compare_references(v2m::vector_type const &ref, std::string_view const &var_ref, std::size_t const var_pos, std::size_t /* out */ &idx)
	{
		char const *var_ref_data(var_ref.data());
		auto const var_ref_len(var_ref.size());
		
		char const *ref_data(ref.data());
		auto const ref_len(ref.size());
		
		if (! (var_pos + var_ref_len <= ref_len))
		{
			idx = 0;
			return false;
		}
		
		auto const var_ref_end(var_ref_data + var_ref_len);
		auto const p(std::mismatch(var_ref_data, var_ref_end, ref_data + var_pos));
		if (var_ref_end != p.first)
		{
			idx = p.first - var_ref_data;
			return false;
		}
		
		return true;
	}
}


namespace vcf2multialign {
	
	void preparation_task::check_ploidy()
	{
		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(vcf_field::ALL);
		
		m_vcf_reader.fill_buffer();
		if (!m_vcf_reader.parse([this](transient_variant const &var) -> bool {
			for (auto const &kv : m_vcf_reader.sample_names())
			{
				auto const sample_no(kv.second);
				auto const &sample(var.sample(sample_no));
				m_ploidy[sample_no] = sample.ploidy();
			}
			
			return false;
		}))
		{
			fail("Unable to read the first variant");
		}
	}
	
	
	void preparation_task::check_ref()
	{
		auto main_queue(dispatch_get_main_queue());

		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(v2m::vcf_field::REF);
		bool found_mismatch(false);
		std::size_t i(0);
		
		bool should_continue(false);
		do {
			m_vcf_reader.fill_buffer();
			should_continue = m_vcf_reader.parse(
				[this, &found_mismatch, &i, main_queue]
				(v2m::transient_variant const &var)
				-> bool
			{
				auto const &var_ref(var.ref());
				auto const var_pos(var.zero_based_pos());
				auto const lineno(var.lineno());
				std::size_t diff_pos{0};
			
				if (!compare_references(*m_reference, var_ref, var_pos, diff_pos))
				{
					if (!found_mismatch)
					{
						found_mismatch = true;
						m_status_logger->log([lineno](){
							std::cerr
								<< "Reference differs from the variant file on line "
								<< lineno
								<< " (and possibly others)."
								<< std::endl;
						});
					}
				
					m_error_logger->log_ref_mismatch(lineno, diff_pos);
				}
				
				return true;
			});
		} while (should_continue);
	}
	
	
	void preparation_task::execute()
	{
		// Update the status in the main queue by calling dispatch_async and synchronize between phases.
		auto main_queue(dispatch_get_main_queue());
		
		// Check the ploidy from the first record.
		m_status_logger->log([](){
			std::cerr << "Checking ploidy…" << std::endl;
		});
		check_ploidy();
		
		// Compare REF to the reference vector.
		// Make the reader count the total number of rows.
		bool did_count(false);
		if (m_should_check_ref)
		{
			m_status_logger->log_message_counting("Comparing the REF column to the reference…");
			check_ref();
			m_record_count = m_vcf_reader.counter_value();
			m_status_logger->finish_logging();
		
			did_count = true;
		}
		
		// List variants that conflict, i.e. overlap but are not nested.
		// Also mark structural variants that cannot be handled.
		auto const msg("Checking overlapping variants…");
		if (did_count)
			m_status_logger->log_message_progress_bar(msg);
		else
			m_status_logger->log_message_counting(msg);
		
		auto const conflict_count(check_overlapping_non_nested_variants(
			m_vcf_reader,
			m_sv_handling_method,
			m_skipped_variants,
			*m_status_logger,
			*m_error_logger
		));
			
		if (!did_count)
			m_record_count = m_vcf_reader.counter_value();
		
		m_status_logger->finish_logging();
		
		dispatch_async_fn(main_queue, [this, did_count, conflict_count]{
			auto const skipped_count(m_skipped_variants.size());
			if (0 == skipped_count)
				std::cerr << "Found no conflicting variants." << std::endl;
			else
			{
				std::cerr << "Found " << conflict_count << " conflicts in total.\n";
				std::cerr << "Number of variants to be skipped: " << m_skipped_variants.size() << std::endl;
			}
		});
		
		// Find subgraphs connected by single edges.
		m_status_logger->log_message_progress_bar("Finding subgraphs connected by single edges…");
		find_subgraph_starting_points(m_vcf_reader, m_skipped_variants, m_subgraph_starting_points);
		m_status_logger->finish_logging();
		
		dispatch_async_fn(main_queue, [this](){
			std::cerr << "Found " << m_subgraph_starting_points.size() << " possible starting points." << std::endl;
			
			// We're done.
			m_delegate->task_did_finish(*this);
		});
		
		// When the block above executes, this will be invalid b.c. m_delegate will deallocate it.
	}
}
