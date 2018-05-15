/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <boost/io/ios_state.hpp>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <libbio/vcf_reader.hh>
#include <map>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>
#include <vcf2multialign/generate_haplotypes.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_handler.hh>


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


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
	
	
	class generate_context
	{
	protected:
		typedef std::map <std::size_t, std::size_t> ploidy_map;

	protected:
		v2m::variant_handler								m_variant_handler;
		lb::dispatch_ptr <dispatch_queue_t>					m_main_queue{};
		lb::dispatch_ptr <dispatch_queue_t>					m_parsing_queue{};
		v2m::vector_type									m_reference;
		lb::vcf_stream_input <lb::file_istream>				m_vcf_input;
		lb::vcf_reader										m_vcf_reader;
		
		std::chrono::time_point <std::chrono::system_clock>	m_start_time{};
		std::chrono::time_point <std::chrono::system_clock>	m_round_start_time{};
		
		v2m::error_logger									m_error_logger;
		
		ploidy_map											m_ploidy;
		v2m::haplotype_map									m_haplotypes;
		v2m::variant_set									m_skipped_variants;

		std::string											m_null_allele_seq;
		lb::vcf_reader::sample_name_map::const_iterator		m_sample_names_it{};
		lb::vcf_reader::sample_name_map::const_iterator		m_sample_names_end{};
		std::size_t											m_chunk_size{0};
		std::size_t											m_current_round{0};
		std::size_t											m_total_rounds{0};
		bool												m_should_overwrite_files{false};
		
	public:
		generate_context(
			lb::dispatch_ptr <dispatch_queue_t> &&main_queue,
			lb::dispatch_ptr <dispatch_queue_t> &&parsing_queue,
			char const *null_allele_seq,
			std::size_t const chunk_size,
			bool const should_overwrite_files
		):
			m_main_queue(std::move(main_queue)),
			m_parsing_queue(std::move(parsing_queue)),
			m_null_allele_seq(null_allele_seq),
			m_chunk_size(chunk_size),
			m_should_overwrite_files(should_overwrite_files)
		{
		}
		
		
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
		
		
		void cleanup() { delete this; }
		
		
		void check_ploidy()
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
				lb::fail("Unable to read the first variant");
			}
		}
		
		
		void check_ref()
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
					auto const var_ref(var.ref());
					auto const var_pos(var.zero_based_pos());
					auto const lineno(var.lineno());
					std::size_t diff_pos{0};
				
					if (!compare_references(m_reference, var_ref, var_pos, diff_pos))
					{
						if (!found_mismatch)
						{
							found_mismatch = true;
							std::cerr << "Reference differs from the variant file on line " << lineno << " (and possibly others)." << std::endl;
						}
					
						m_error_logger.log_ref_mismatch(lineno, diff_pos);
					}
					
					++i;
					if (0 == i % 100000)
						std::cerr << "Handled " << i << " variants…" << std::endl;
					
					return true;
				});
			} while (should_continue);
		}
		
		
		void update_haplotypes(char const *out_reference_fname)
		{
			m_haplotypes.clear();
			
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
					lb::open_file_for_writing(fname.c_str(), haplotype_vec[i - 1].output_stream, m_should_overwrite_files);
				}

				++m_sample_names_it;
				++i;
				
				if (i == m_chunk_size)
					break;
			}
			
			// Check if reference output was requested.
			if (out_reference_fname)
			{
				lb::always_assert(
					m_haplotypes.cend() == m_haplotypes.find(v2m::REF_SAMPLE_NUMBER),
					"REF_SAMPLE_NUMBER already in use"
				);
				
				auto it(m_haplotypes.emplace(
					std::piecewise_construct,
					std::forward_as_tuple(0),
					std::forward_as_tuple(1)
				).first);
				
				auto &haplotype_vec(it->second);
				lb::open_file_for_writing(out_reference_fname, haplotype_vec[0].output_stream, m_should_overwrite_files);
			}
		}

		
		// Handle m_chunk_size samples.
		void generate_sequences(char const *out_reference_fname = nullptr)
		{
			// Open files for the samples. If no files were opened, exit.
			update_haplotypes(out_reference_fname);
			if (0 == m_haplotypes.size())
			{
				{
					// Save the stream state.
					boost::io::ios_flags_saver ifs(std::cerr);
				
					// Change FP notation.
					std::cerr << std::fixed;

					auto const end_time(std::chrono::system_clock::now());
					std::chrono::duration <double> elapsed_seconds(end_time - m_start_time);
					std::cerr << "Sequence generation took " << (elapsed_seconds.count() / 60.0) << " minutes in total." << std::endl;
				}

				// After calling cleanup *this is no longer valid.
				//std::cerr << "Calling cleanup" << std::endl;
				cleanup();
				
				//std::cerr << "Calling exit" << std::endl;
				exit(EXIT_SUCCESS);
			}
			
			++m_current_round;
			std::cerr << "Round " << m_current_round << '/' << m_total_rounds << std::endl;
			m_round_start_time = std::chrono::system_clock::now();
			auto const start_time(std::chrono::system_clock::to_time_t(m_round_start_time));
			std::cerr << "Starting on " << std::ctime(&start_time) << std::flush;
			m_variant_handler.process_variants(m_haplotypes);
		}
		
		
		void finish_round()
		{
			// Save the stream state.
			boost::io::ios_flags_saver ifs(std::cerr);

			// Change FP notation.
			std::cerr << std::fixed;
			
			auto const end_time(std::chrono::system_clock::now());
			std::chrono::duration <double> elapsed_seconds(end_time - m_round_start_time);
			std::cerr << "Finished in " << (elapsed_seconds.count() / 60.0) << " minutes." << std::endl;
		}

		
		void load_and_generate(
			char const *reference_fname,
			char const *variants_fname,
			char const *out_reference_fname,
			char const *report_fname,
			v2m::sv_handling const sv_handling_method,
			bool const should_check_ref
		)
		{
			// Open the files.
			std::cerr << "Opening files…" << std::endl;
			{
				lb::file_istream ref_fasta_stream;
				
				lb::open_file_for_reading(reference_fname, ref_fasta_stream);
				lb::open_file_for_reading(variants_fname, m_vcf_input.input_stream());
				
				if (report_fname)
				{
					lb::open_file_for_writing(report_fname, m_error_logger.output_stream(), m_should_overwrite_files);
					m_error_logger.write_header();
				}
				
				m_vcf_reader.set_input(m_vcf_input);
				m_vcf_reader.read_header();
				
				auto const &sample_names(m_vcf_reader.sample_names());
				m_sample_names_it = sample_names.cbegin();
				m_sample_names_end = sample_names.cend();

				auto const sample_count(sample_names.size());
				m_total_rounds = std::ceil(1.0 * sample_count / m_chunk_size);
				
				// Read the reference file and place its contents into reference.
				v2m::read_single_fasta_seq(ref_fasta_stream, m_reference);
			}
			
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
				auto const conflict_count(v2m::check_overlapping_non_nested_variants(
					m_vcf_reader,
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
			
			// Replace the placeholder variant_handler.
			{
				v2m::variant_handler temp(
					m_main_queue,
					m_parsing_queue,
					m_vcf_reader,
					m_reference,
					sv_handling_method,
					m_skipped_variants,
					m_null_allele_seq,
					m_error_logger,
					[this](){ finish_round(); generate_sequences(); }
				);
				
				m_variant_handler = std::move(temp);
				m_variant_handler.get_variant_buffer().set_delegate(m_variant_handler);
			}
			
			std::cerr << "Generating haplotype sequences…" << std::endl;
			m_start_time = std::chrono::system_clock::now();
			generate_sequences(out_reference_fname);
		}
	};
}


namespace vcf2multialign {
	
	void generate_haplotypes(
		char const *reference_fname,
		char const *variants_fname,
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
		
		// generate_context needs to be allocated on the heap because later dispatch_main is called.
		// The class deallocates itself in cleanup().
		generate_context *ctx(new generate_context(
			std::move(main_queue),
			std::move(parsing_queue),
			null_allele_seq,
			chunk_size,
			should_overwrite_files
		));
		
		ctx->load_and_generate(
			reference_fname,
			variants_fname,
			out_reference_fname,
			report_fname,
			sv_handling_method,
			should_check_ref
		);
	}
}
