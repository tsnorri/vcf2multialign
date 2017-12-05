/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/io/ios_state.hpp>
#include <vcf2multialign/tasks/all_haplotypes_task.hh>


namespace vcf2multialign {
	
	auto all_haplotypes_task::create_haplotype(
		std::size_t const sample_no,
		std::size_t const ploidy
	) -> haplotype_map_type::iterator
	{
		return m_haplotypes.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(sample_no),
			std::forward_as_tuple(ploidy)
		).first;
	}
	
	
	auto all_haplotypes_task::find_or_create_haplotype(
		std::size_t const sample_no,
		std::size_t const ploidy
	) -> haplotype_map_type::iterator
	{
		// Since file_ostream is not movable, check first if the vector has already been created.
		// If not, create it with the exact size to avoid resizing and thus moving later.
		auto it(m_haplotypes.find(sample_no));
		if (m_haplotypes.end() == it)
			it = create_haplotype(sample_no, ploidy);
		
		return it;
	}
	
	
	void all_haplotypes_task::update_haplotypes_with_ref()
	{
		assert(m_out_reference_fname->operator bool());
		
		// Check if reference output was requested.
		always_assert(
			m_haplotypes.cend() == m_haplotypes.find(REF_SAMPLE_NUMBER),
			"REF_SAMPLE_NUMBER already in use"
		);
		
		auto it(create_haplotype(REF_SAMPLE_NUMBER, 1));
		auto &haplotype_vec(it->second);
		open_file_for_writing(
			(*m_out_reference_fname)->c_str(),
			haplotype_vec[0].output_stream,
			m_should_overwrite_files
		);
	}
	
	
	void all_haplotypes_task::update_haplotypes(bool const output_reference)
	{
		std::size_t i(0);
		while (m_sample_names_it != m_sample_names_end)
		{
			auto const &sample_name(m_sample_names_it->first);
			auto const sample_no(m_sample_names_it->second);
			auto const current_ploidy(m_ploidy->find(sample_no)->second);
		
			auto it(find_or_create_haplotype(sample_no, current_ploidy));
			auto &haplotype_vec(it->second);
		
			for (std::size_t j(1); j <= current_ploidy; ++j)
			{
				auto const fname(boost::str(boost::format("%s-%u") % sample_name % j));
				open_file_for_writing(fname.c_str(), haplotype_vec[j - 1].output_stream, m_should_overwrite_files);
			}
	
			++m_sample_names_it;
			++i;
		
			if (i == m_chunk_size - (output_reference ? 1 : 0))
				break;
		}
		
		if (output_reference)
			update_haplotypes_with_ref();
	}
	
	
	void all_haplotypes_task::generate_sequences(bool const output_reference)
	{
		// Open files for the samples. If no files were opened, exit.
		m_haplotypes.clear();
		update_haplotypes(output_reference);
	
		auto const haplotype_count(m_haplotypes.size());
		if (0 == haplotype_count)
		{
			auto const start_time(m_start_time);
			m_logger->status_logger.log([start_time](){
				// Save the stream state.
				boost::io::ios_flags_saver ifs(std::cerr);
		
				// Change FP notation.
				std::cerr << std::fixed;
			
				auto const end_time(std::chrono::system_clock::now());
				std::chrono::duration <double> elapsed_seconds(end_time - start_time);
				std::cerr << "Sequence generation took " << (elapsed_seconds.count() / 60.0) << " minutes in total." << std::endl;
			});
			
			auto main_queue(dispatch_get_main_queue());
			dispatch_async_fn(main_queue, [this](){
				m_delegate->task_did_finish(*this);
			});
			
			return;
		}
	
		++m_current_round;
		m_round_start_time = std::chrono::system_clock::now();
	
		auto const round_start_time(m_round_start_time);
		auto const current_round(m_current_round);
		auto const total_rounds(m_total_rounds);
		m_logger->status_logger.log([round_start_time, current_round, total_rounds](){
			auto const start_time(std::chrono::system_clock::to_time_t(round_start_time));
			std::cerr << "Round " << current_round << '/' << total_rounds << std::endl;
			std::cerr << "Starting on " << std::ctime(&start_time) << std::flush;
		});
		m_logger->status_logger.log_message_progress_bar("Generating haplotypes…");
		
		m_sequence_writer.prepare(m_haplotypes);
		variant_handler().process_variants();
	}
	
	
	void all_haplotypes_task::handle_variant(variant &var)
	{
		variant_stats::handle_variant(var);
		m_sequence_writer.handle_variant(var);
	}
	
	
	void all_haplotypes_task::finish()
	{
		m_sequence_writer.finish();
		m_logger->status_logger.finish_logging();
		
		auto const round_start_time(m_round_start_time);
		m_logger->status_logger.log([round_start_time](){
			// Save the stream state.
			boost::io::ios_flags_saver ifs(std::cerr);
	
			// Change FP notation.
			std::cerr << std::fixed;
		
			auto const end_time(std::chrono::system_clock::now());
			std::chrono::duration <double> elapsed_seconds(end_time - round_start_time);
			std::cerr << "Finished in " << (elapsed_seconds.count() / 60.0) << " minutes." << std::endl;
		});
		
		generate_sequences();
	}
	
	
	void all_haplotypes_task::execute()
	{
		auto &vr(vcf_reader());
		
		// Prepare sample names for enumeration and calculate rounds.
		auto const &sample_names(vr.sample_names());
		auto const sample_count(sample_names.size());
		m_total_rounds = std::ceil(1.0 * sample_count / m_chunk_size);
		
		m_sample_names_it = sample_names.cbegin();
		m_sample_names_end = sample_names.cend();
		
		m_logger->status_logger.log([](){
			std::cerr << "Generating haplotype sequences…" << std::endl;
		});

		vr.set_parsed_fields(vcf_field::ALL);

		m_start_time = std::chrono::system_clock::now();
		generate_sequences(m_out_reference_fname->operator bool());
	}
}
