/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/io/ios_state.hpp>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/generate_context_base.hh>


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	bool compare_references(
		v2m::vector_type const &ref,
		std::string_view const &var_ref,
		std::size_t const var_pos,
		std::size_t /* out */ &idx
	)
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
	
	void generate_context_base::open_files(
		char const *reference_fname,
		char const *ref_seq_name,
		char const *variants_fname,
		char const *report_fname
	)
	{
		lb::mmap_handle <char> ref_handle;
		ref_handle.open(reference_fname);
		
		lb::open_file_for_reading(variants_fname, m_vcf_input.input_stream());
		
		if (report_fname)
		{
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(m_should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			
			lb::open_file_for_writing(report_fname, m_error_logger.output_stream(), mode);
			m_error_logger.write_header();
		}
		
		m_vcf_reader.set_input(m_vcf_input);
		m_vcf_reader.read_header();
		
		// Read the reference file and place its contents into reference.
		read_single_fasta_seq(ref_handle, m_reference, ref_seq_name);
	}
	
		
	void generate_context_base::check_ploidy()
	{
		size_t i(0);
		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(lb::vcf_field::ALL);
		
		m_vcf_reader.fill_buffer();
		if (!m_vcf_reader.parse([this](lb::transient_variant const &var) -> bool {
			for (auto const &kv : m_vcf_reader.sample_names())
			{
				auto const sample_no(kv.second);
				auto const &sample(var.sample(sample_no));
				m_ploidy[sample_no] = sample.ploidy();
			}
			
			return false;
		}))
		{
			libbio_fail("Unable to read the first variant");
		}
	}
	
	
	void generate_context_base::check_ref()
	{
		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(lb::vcf_field::REF);
		bool found_mismatch(false);
		std::size_t i(0);
		
		bool should_continue(false);
		do {
			m_vcf_reader.fill_buffer();
			should_continue = m_vcf_reader.parse(
				[this, &found_mismatch, &i]
				(lb::transient_variant const &var)
				-> bool
			{
				if (0 == m_chromosome_name.size() || var.chrom_id() == m_chromosome_name)
				{
					auto const var_ref(var.ref());
					auto const var_pos(var.zero_based_pos());
					auto const lineno(var.lineno());
					std::size_t diff_pos{0};
			
					if (!compare_references(m_reference, var_ref, var_pos, diff_pos))
					{
						if (!found_mismatch)
						{
							found_mismatch = true;
							std::cerr
								<< "Reference differs from the variant file on line "
								<< lineno
								<< " (and possibly others)."
								<< std::endl;
						}
				
						m_error_logger.log_ref_mismatch(lineno, diff_pos);
					}
				}
				
				++i;
				if (0 == i % 100000)
					std::cerr << "Handled " << i << " variants…" << std::endl;
				
				return true;
			});
		} while (should_continue);
	}
	
	
	void generate_context_base::finish()
	{
		// Save the stream state.
		boost::io::ios_flags_saver ifs(std::cerr);
		
		// Change FP notation.
		std::cerr << std::fixed;
		
		auto const end_time(std::chrono::system_clock::now());
		std::chrono::duration <double> elapsed_seconds(end_time - m_start_time);
		std::cerr << "Sequence generation took " << (elapsed_seconds.count() / 60.0) << " minutes in total." << std::endl;
		
		// After calling cleanup *this is no longer valid.
		//std::cerr << "Calling cleanup" << std::endl;
		cleanup();
		
		//std::cerr << "Calling exit" << std::endl;
		exit(EXIT_SUCCESS);
	}
	
	
	void generate_context_base::load_and_generate(
		sv_handling const sv_handling_method,
		bool const should_check_ref
	)
	{
		// Check ploidy from the first record.
		std::cerr << "Checking ploidy…" << std::endl;
		check_ploidy();
		
		// Compare REF to the reference vector.
		if (should_check_ref)
		{
			std::cerr << "Comparing the REF column to the reference…" << std::endl;
			check_ref();
		}
		
		// List variants that conflict, i.e. overlap but are not nested.
		// Also mark structural variants that cannot be handled.
		{
			std::cerr << "Checking overlapping variants…" << std::endl;
			auto const conflict_count(check_overlapping_non_nested_variants(
				m_vcf_reader,
				m_chromosome_name,
				sv_handling_method,
				m_skipped_variants,
				m_error_logger
			));

			auto const skipped_count(m_skipped_variants.size());
			if (0 == skipped_count)
				std::cerr << "Found no conflicting variants." << std::endl;
			else
			{
				std::cerr << "Found " << conflict_count << " conflicts in total." << std::endl;
				std::cerr << "Number of variants to be skipped: " << m_skipped_variants.size() << std::endl;
			}
		}
	}
}
