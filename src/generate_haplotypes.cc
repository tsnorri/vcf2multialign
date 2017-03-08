/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <map>
#include <vcf2multialign/generate_haplotypes.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_handler.hh>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>

namespace ios	= boost::iostreams;
namespace v2m	= vcf2multialign;


namespace {

	void handle_file_error(char const *fname)
	{
		char const *errmsg(strerror(errno));
		std::cerr << "Got an error while trying to open '" << fname << "': " << errmsg << std::endl;
		exit(EXIT_FAILURE);
	}


	void open_file_for_reading(char const *fname, v2m::file_istream &stream)
	{
		int fd(open(fname, O_RDONLY));
		if (-1 == fd)
			handle_file_error(fname);

		ios::file_descriptor_source source(fd, ios::close_handle);
		stream.open(source);
		stream.exceptions(std::istream::badbit);
	}
	
	
	void open_file_for_writing(char const *fname, v2m::file_ostream &stream, bool const should_overwrite)
	{
		int fd(0);
		if (should_overwrite)
			fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		else
			fd = open(fname, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		
		if (-1 == fd)
			handle_file_error(fname);
		
		ios::file_descriptor_sink sink(fd, ios::close_handle);
		stream.open(sink);
	}
	
	
	class generate_context
	{
	protected:
		typedef std::map <std::size_t, std::size_t> ploidy_map;

	protected:
		v2m::variant_handler								m_variant_handler;
		dispatch_queue_t									m_main_queue{};
		dispatch_queue_t									m_parsing_queue{};
		v2m::vector_type									m_reference;
		v2m::file_istream									m_vcf_stream;
		v2m::vcf_reader										m_vcf_reader;
		
		v2m::error_logger									m_error_logger;
		
		ploidy_map											m_ploidy;
		v2m::haplotype_map									m_haplotypes;
		v2m::variant_set									m_skipped_variants;

		std::string											m_null_allele_seq;
		v2m::vcf_reader::sample_name_map::const_iterator	m_sample_names_it{};
		v2m::vcf_reader::sample_name_map::const_iterator	m_sample_names_end{};
		std::size_t											m_chunk_size{0};
		std::size_t											m_current_round{0};
		std::size_t											m_total_rounds{0};
		bool												m_should_overwrite_files{false};
		
	public:
		generate_context(
			dispatch_queue_t main_queue,
			dispatch_queue_t parsing_queue,
			char const *null_allele_seq,
			std::size_t const chunk_size,
			bool const should_overwrite_files
		):
			m_main_queue(main_queue),
			m_parsing_queue(parsing_queue),
			m_null_allele_seq(null_allele_seq),
			m_chunk_size(chunk_size),
			m_should_overwrite_files(should_overwrite_files)
		{
			dispatch_retain(m_main_queue);
			dispatch_retain(m_parsing_queue);
		}
		
		
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
		
		
		~generate_context()
		{
			dispatch_release(m_main_queue);
			dispatch_release(m_parsing_queue);
		}
		
		
		void cleanup() { delete this; }
		
		
		void check_ploidy()
		{
			size_t i(0);
			m_vcf_reader.reset();
			m_vcf_reader.set_parsed_fields(v2m::vcf_field::ALL);
			v2m::variant var(m_vcf_reader.sample_count());
			var.add_format_field("GT");

			if (!m_vcf_reader.get_next_variant(var))
				throw std::runtime_error("Unable to read the first variant");
			
			var.prepare_samples();
			v2m::variant::sample_field_vector sample_fields;
			v2m::variant::genotype_vector gtv;
			bool is_phased(false);
			
			for (auto const &kv : m_vcf_reader.sample_names())
			{
				auto const sample_no(kv.second);
				var.parse_sample(sample_no, sample_fields);
				var.get_genotype(sample_fields, gtv, is_phased);
				m_ploidy[sample_no] = gtv.size();
			}
		}
		
		
		void update_haplotypes()
		{
			m_haplotypes.clear();
			
			size_t i(1);
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
					open_file_for_writing(fname.c_str(), haplotype_vec[i - 1].output_stream, m_should_overwrite_files);
				}

				++m_sample_names_it;
				++i;
				
				if (i == m_chunk_size)
					break;
			}
		}

		
		// Handle m_chunk_size samples.
		void generate_sequences()
		{
			// Open files for the samples. If no files were opened, exit.
			update_haplotypes();
			if (0 == m_haplotypes.size())
			{
				// After calling cleanup *this is no longer valid.
				//std::cerr << "Calling cleanup" << std::endl;
				cleanup();
				
				//std::cerr << "Calling exit" << std::endl;
				exit(EXIT_SUCCESS);
			}
			
			++m_current_round;
			std::cerr << "Round " << m_current_round << '/' << m_total_rounds << std::endl;
			m_variant_handler.process_variants(m_haplotypes);
		}

		
		void load_and_generate(
			char const *reference_fname,
			char const *variants_fname,
			char const *report_fname
		)
		{
			// Open the files.
			{
				v2m::file_istream ref_fasta_stream;
				
				open_file_for_reading(reference_fname, ref_fasta_stream);
				open_file_for_reading(variants_fname, m_vcf_stream);
				
				if (report_fname)
				{
					open_file_for_writing(report_fname, m_error_logger.output_stream(), m_should_overwrite_files);
					m_error_logger.write_header();
				}
				
				m_vcf_reader.set_stream(m_vcf_stream);
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
			check_ploidy();
			
			// List variants that conflict, i.e. overlap but are not nested.
			{
				auto const conflict_count(v2m::check_overlapping_non_nested_variants(m_vcf_reader, m_skipped_variants, m_error_logger));
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
					m_skipped_variants,
					m_null_allele_seq,
					m_error_logger,
					[this](){ generate_sequences(); }
				);
				
				m_variant_handler = std::move(temp);
			}
			
			std::cerr << "Generating haplotype sequences…" << std::endl;
			generate_sequences();
		}
	};
}


namespace vcf2multialign {
	
	void generate_haplotypes(
		char const *reference_fname,
		char const *variants_fname,
		char const *report_fname,
		char const *null_allele_seq,
		std::size_t const chunk_size,
		bool const should_overwrite_files
	)
	{
		dispatch_queue_t main_queue(dispatch_get_main_queue());
		dispatch_queue_t parsing_queue(dispatch_queue_create("fi.iki.tsnorri.vcf2multialign.parsing_queue", DISPATCH_QUEUE_CONCURRENT));
		
		// generate_context needs to be allocated on the heap because later dispatch_main is called.
		// The class deallocates itself in cleanup().
		generate_context *ctx(new generate_context(
			main_queue,
			parsing_queue,
			null_allele_seq,
			chunk_size,
			should_overwrite_files
		));
		
		ctx->load_and_generate(
			reference_fname,
			variants_fname,
			report_fname
		);
		
		// Calls pthread_exit.
		dispatch_main();
	}
}
