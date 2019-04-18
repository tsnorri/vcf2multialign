/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <boost/io/ios_state.hpp>
#include <vcf2multialign/sequences/generate_sequences_context.hh>


namespace lb = libbio;


namespace vcf2multialign {
	
	void generate_sequences_context::update_haplotypes(char const *out_reference_fname)
	{
		m_haplotypes.clear();
		
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(m_should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		
		size_t i(0);
		while (m_sample_names_it != m_sample_names_end)
		{
			auto const &sample_name(m_sample_names_it->first);
			auto const sample_no(m_sample_names_it->second);
			auto const current_ploidy(m_ploidy.find(sample_no)->second);
			
			// Since file_ostream is not movable, check first if the vector has already been created.
			// If not, create it with the exact size to avoid resizing and thus moving later.
			auto it(m_haplotypes.find(sample_no));
			if (m_haplotypes.end() == it)
			{
				it = m_haplotypes.emplace(
					std::piecewise_construct,
					std::forward_as_tuple(sample_no),
					std::forward_as_tuple(current_ploidy)
				).first;
			}
			
			auto &haplotype_vec(it->second);
			
			for (size_t i(1); i <= current_ploidy; ++i)
			{
				auto const fname(boost::str(boost::format("%s-%u") % sample_name % i));
				lb::open_file_for_writing(fname.c_str(), haplotype_vec[i - 1].output_stream, mode);
			}
	
			++m_sample_names_it;
			++i;
			
			if (i == m_chunk_size)
				break;
		}
		
		// Check if reference output was requested.
		if (out_reference_fname)
		{
			libbio_always_assert_msg(
				m_haplotypes.cend() == m_haplotypes.find(REF_SAMPLE_NUMBER),
				"REF_SAMPLE_NUMBER already in use"
			);
			
			auto it(m_haplotypes.emplace(
				std::piecewise_construct,
				std::forward_as_tuple(REF_SAMPLE_NUMBER),
				std::forward_as_tuple(1)
			).first);
			
			auto &haplotype_vec(it->second);
			lb::open_file_for_writing(out_reference_fname, haplotype_vec[0].output_stream, mode);
		}
	}
	
	
	// Handle m_chunk_size samples.
	void generate_sequences_context::generate_sequences(char const *out_reference_fname)
	{
		// Open files for the samples. If no files were opened, exit.
		update_haplotypes(out_reference_fname);
		if (0 == m_haplotypes.size())
			finish();
		
		++m_current_round;
		std::cerr << "Round " << m_current_round << '/' << m_total_rounds << std::endl;
		m_round_start_time = std::chrono::system_clock::now();
		auto const start_time(std::chrono::system_clock::to_time_t(m_round_start_time));
		std::cerr << "Starting on " << std::ctime(&start_time) << std::flush;
		m_variant_handler.process_variants(m_haplotypes);
	}
	
	
	void generate_sequences_context::finish_round()
	{
		// Save the stream state.
		boost::io::ios_flags_saver ifs(std::cerr);
	
		// Change FP notation.
		std::cerr << std::fixed;
		
		auto const end_time(std::chrono::system_clock::now());
		std::chrono::duration <double> elapsed_seconds(end_time - m_round_start_time);
		std::cerr << "Finished in " << (elapsed_seconds.count() / 60.0) << " minutes." << std::endl;
	}
	
	
	void generate_sequences_context::variant_handler_did_finish(variant_handler_type &handler)
	{
		finish_round();
		generate_sequences();
	}
	
	
	void generate_sequences_context::load_and_generate(
		char const *reference_fname,
		char const *ref_seq_name,
		char const *variants_fname,
		char const *out_reference_fname,
		char const *report_fname,
		bool const should_check_ref
	)
	{
		// Open the files.
		std::cerr << "Opening the files…" << std::endl;
		open_files(reference_fname, ref_seq_name, variants_fname, report_fname);
		
		{
			auto const &sample_names(m_vcf_reader.sample_names());
			m_sample_names_it = sample_names.cbegin();
			m_sample_names_end = sample_names.cend();

			auto const sample_count(sample_names.size());
			libbio_always_assert(sample_count);
			m_total_rounds = std::ceil(double(sample_count) / m_chunk_size);
		}
		
		generate_context_base::load_and_generate(should_check_ref);
			
		// Replace the placeholder variant_handler_type.
		{
			variant_handler_type temp(
				m_main_queue,
				m_parsing_queue,
				m_vcf_reader,
				m_reference,
				m_chromosome_name,
				m_skipped_variants,
				m_null_allele_seq,
				m_error_logger,
				*this
			);
			
			m_variant_handler = std::move(temp);
			m_variant_handler.get_variant_buffer().set_delegate(m_variant_handler);
		}
		
		std::cerr << "Generating haplotype sequences…" << std::endl;
		m_start_time = std::chrono::system_clock::now();
		generate_sequences(out_reference_fname);
	}
}
