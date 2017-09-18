/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/range/combine.hpp>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <map>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/generate_haplotypes.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/sample_reducer.hh>
#include <vcf2multialign/sequence_writer.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_handler.hh>

namespace ios	= boost::iostreams;
namespace v2m	= vcf2multialign;


namespace {
	class generate_context;
	class genotype_handling_delegate;
	class all_genotypes_handling_delegate;
	class compressed_genotypes_handling_delegate;
	
	void handle_file_error(char const *fname);
	void open_file_for_reading(char const *fname, v2m::file_istream &stream);
	void open_file_for_writing(char const *fname, v2m::file_ostream &stream, bool const should_overwrite);
	bool compare_references(v2m::vector_type const &ref, std::string_view const &var_ref, std::size_t const var_pos, std::size_t /* out */ &idx);
	
	v2m::haplotype_map::iterator create_haplotype(
		v2m::haplotype_map &haplotypes,
		std::size_t const sample_no,
		std::size_t const ploidy
	);
	
	v2m::haplotype_map::iterator find_or_create_haplotype(
		v2m::haplotype_map &haplotypes,
		std::size_t const sample_no,
		std::size_t const ploidy
	);
		
		
	class generate_context
	{
	protected:
		typedef std::map <std::size_t, std::size_t> ploidy_map;
	
	protected:
		v2m::variant_handler								m_variant_handler;
	
		v2m::vector_type									m_reference;
		v2m::file_istream									m_vcf_stream;
		v2m::vcf_reader										m_vcf_reader;
	
		std::unique_ptr <v2m::variant_handler_delegate>		m_variant_handler_delegate{};
		std::unique_ptr <genotype_handling_delegate>		m_genotype_delegate{};
	
		std::chrono::time_point <std::chrono::system_clock>	m_start_time{};
		std::chrono::time_point <std::chrono::system_clock>	m_round_start_time{};
	
		v2m::error_logger									m_error_logger;
	
		ploidy_map											m_ploidy;
		v2m::haplotype_map									m_haplotypes;
		v2m::variant_set									m_skipped_variants;
	
		boost::optional <std::string>						m_out_reference_fname;
		std::string											m_null_allele_seq;
		v2m::sv_handling									m_sv_handling_method;
		std::size_t											m_chunk_size{0};
		std::size_t											m_current_round{0};
		std::size_t											m_total_rounds{0};
		bool												m_should_overwrite_files{false};
	
	public:
		generate_context(
			v2m::dispatch_ptr <dispatch_queue_t> &&main_queue,
			v2m::dispatch_ptr <dispatch_queue_t> &&parsing_queue,
			char const *out_reference_fname,
			char const *null_allele_seq,
			v2m::sv_handling const sv_handling_method,
			std::size_t const chunk_size,
			std::size_t const variant_padding,
			bool const should_overwrite_files,
			bool const should_reduce_samples,
			bool const allow_switch_to_ref
		):
			m_variant_handler(
				std::move(main_queue),
				std::move(parsing_queue),
				m_vcf_reader,
				m_reference,
				sv_handling_method,
				m_skipped_variants,
				m_error_logger
			),
			m_null_allele_seq(null_allele_seq),
			m_sv_handling_method(sv_handling_method),
			m_chunk_size(chunk_size),
			m_should_overwrite_files(should_overwrite_files)
		{
			finish_init(
				out_reference_fname,
				should_reduce_samples,
				variant_padding,
				allow_switch_to_ref
			);
		}
	
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
	
		v2m::vcf_reader &vcf_reader()						{ return m_vcf_reader; }
		v2m::error_logger &error_logger()					{ return m_error_logger; }
		v2m::variant_handler &variant_handler()				{ return m_variant_handler; }
		v2m::haplotype_map &haplotypes()					{ return m_haplotypes; }

		v2m::haplotype_map const &haplotypes() const		{ return m_haplotypes; }
		v2m::variant_handler const &variant_handler() const	{ return m_variant_handler; }
		ploidy_map const &ploidy() const					{ return m_ploidy; }
		v2m::variant_set const &skipped_variants() const	{ return m_skipped_variants; }
		std::size_t chunk_size() const						{ return m_chunk_size; }
		std::string const &out_reference_fname() const		{ return m_out_reference_fname.value(); }
		std::string const &null_allele_seq() const			{ return m_null_allele_seq; }
		v2m::vector_type const &reference() const			{ return m_reference; }
		bool should_overwrite_files() const					{ return m_should_overwrite_files; }
		bool has_out_reference_fname() const				{ return m_out_reference_fname.operator bool(); }
		
		void set_variant_handler_delegate(std::unique_ptr <v2m::variant_handler_delegate> &&delegate);
		
		void cleanup() { delete this; }
		void load_and_generate(
			char const *reference_fname,
			char const *variants_fname,
			char const *report_fname,
			bool const should_check_ref
		);
			
		void prepare_sample_names_and_generate_sequences();
		void generate_sequences(bool const output_reference = false);
		void finish_round();
		void process_variants() { m_variant_handler.process_variants(); }
	
	protected:
		void finish_init(
			char const *out_reference_fname,
			bool const should_reduce_samples,
			std::size_t const variant_padding,
			bool const allow_switch_to_ref
		);
		void check_ploidy();
		void check_ref();
	};
	
	
	class genotype_handling_delegate
	{
	protected:
		generate_context *m_generate_context{};
		
	public:
		genotype_handling_delegate() = default;
		virtual ~genotype_handling_delegate() {}
		
		void set_generate_context(generate_context &ctx) { m_generate_context = &ctx; }
		void update_haplotypes_with_ref(v2m::haplotype_map &haplotypes);

		virtual void process_first_phase() = 0;
		virtual void prepare_sample_names() = 0;
		virtual std::size_t sample_count() const = 0;
		virtual void update_haplotypes(v2m::haplotype_map &haplotypes, bool const output_reference) = 0;
	};
	
	
	class all_genotypes_handling_delegate : public genotype_handling_delegate
	{
	protected:
		v2m::vcf_reader::sample_name_map::const_iterator	m_sample_names_it{};
		v2m::vcf_reader::sample_name_map::const_iterator	m_sample_names_end{};
		std::size_t											m_sample_count{};
		
	public:
		void process_first_phase() override;
		void prepare_sample_names() override;
		std::size_t sample_count() const override { return m_sample_count; }
		void update_haplotypes(v2m::haplotype_map &haplotypes, bool const output_reference) override;
	};
	
	
	class compressed_genotypes_handling_delegate : public genotype_handling_delegate
	{
	protected:
		v2m::range_map	m_compressed_ranges{};
		std::size_t		m_sample_idx{0};
		std::size_t		m_variant_padding{0};
		bool			m_output_reference{};
		bool			m_allow_switch_to_ref{};
		
	public:
		compressed_genotypes_handling_delegate(
			bool output_reference,
			std::size_t const variant_padding,
			bool const allow_switch_to_ref
		):
			m_variant_padding(variant_padding),
			m_output_reference(output_reference),
			m_allow_switch_to_ref(allow_switch_to_ref)
		{
		}
		void process_first_phase() override;
		void prepare_sample_names() override {}
		std::size_t sample_count() const override { return m_compressed_ranges.size(); }
		void update_haplotypes(v2m::haplotype_map &haplotypes, bool const output_reference) override;
	};
	
	
	struct vh_generate_context_helper
	{
		virtual class generate_context &generate_context() = 0;
		virtual class generate_context const &generate_context() const = 0;
	};
	
	
	template <bool t_logging>
	class vh_stats : public virtual v2m::variant_processor_delegate
	{
	public:
		virtual void assigned_alt_to_sequence(std::size_t const alt_idx) override {}
		virtual void found_overlapping_alt(
			std::size_t const lineno,
			uint8_t const alt_idx,
			std::size_t const sample_no,
			uint8_t const chr_idx
		) override {}
		virtual void handled_alt(std::size_t const alt_idx) override {}
		virtual void handled_haplotypes(v2m::variant &var) override {}
	};
	
	
	template <>
	class vh_stats <true> : public virtual v2m::variant_processor_delegate, public virtual vh_generate_context_helper
	{
	protected:
		v2m::variant_set						m_overlapping_alts{};
		std::vector <v2m::skipped_sample>		m_skipped_samples;		// In current variant.
		std::map <uint8_t, v2m::sample_count>	m_counts_by_alt;		// In current variant.
		v2m::sample_count						m_non_ref_totals;		// In current variant.
		
	public:
		virtual void assigned_alt_to_sequence(std::size_t const alt_idx) override;
		virtual void found_overlapping_alt(
			std::size_t const lineno,
			uint8_t const alt_idx,
			std::size_t const sample_no,
			uint8_t const chr_idx
		) override;
		virtual void handled_alt(std::size_t const alt_idx) override;
		virtual void handled_haplotypes(v2m::variant &var) override;

		void handle_variant(v2m::variant &var);
	};
	
	
	class vh_sequence_writer :
		public virtual vh_generate_context_helper,
		public virtual v2m::variant_handler_delegate
	{
	protected:
		v2m::sequence_writer		m_sequence_writer;
		
	public:
		virtual ~vh_sequence_writer() {}
		vh_sequence_writer(
			v2m::sequence_writer_delegate &delegate,
			v2m::vector_type const &reference,
			std::string const &null_allele
		):
			m_sequence_writer(reference, null_allele)
		{
			m_sequence_writer.set_delegate(delegate);
		}
		
		virtual void finish();
	};
	
	
	class vh_delegate :
		public v2m::sample_reducer_delegate,
		public v2m::sequence_writer_delegate,
		public virtual v2m::variant_handler_delegate,
		public virtual vh_generate_context_helper
	{
	protected:
		class generate_context			*m_ctx{};
		
	public:
		virtual ~vh_delegate() {}
		
		vh_delegate(class generate_context &ctx):
			m_ctx(&ctx)
		{
		}
		
		virtual class generate_context &generate_context() override { return *m_ctx; }
		virtual class generate_context const &generate_context() const override { return *m_ctx; }
		
		virtual v2m::vcf_reader::sample_name_map const &sample_names() const override;
		virtual bool is_valid_alt(std::size_t const alt_idx) const override;
		virtual std::set <std::size_t> const &valid_alts(v2m::variant &var) const override;
		
		virtual void enumerate_genotype(
			v2m::variant &var,
			std::size_t const sample_no,
			std::function <void(uint8_t, std::size_t, bool)> const &cb
		) override;
	};
	
	
	class compress_vh_delegate : public vh_stats <true>, public vh_delegate
	{
	protected:
		v2m::sample_reducer		m_sample_reducer;
		
	public:
		compress_vh_delegate(
			class generate_context	&ctx,
			v2m::range_map			&compressed_ranges,
			std::size_t const		padding_amt,
			bool const				output_ref,
			bool const				allow_switch_to_ref
		):
			vh_delegate(ctx),
			m_sample_reducer(compressed_ranges, padding_amt, output_ref, allow_switch_to_ref)
		{
			m_sample_reducer.set_delegate(*this);
		}
		
		virtual void prepare(v2m::vcf_reader &reader) override;
		virtual void handle_variant(v2m::variant &var) override;
		virtual void finish() override;
	};
	
	
	class read_compressed_vh_delegate : public vh_stats <false>, public vh_sequence_writer, public vh_delegate
	{
	protected:
		typedef v2m::range_map::value_type::const_iterator		variant_range_iterator;
		typedef boost::iterator_range <variant_range_iterator>	range_el_iterator;
		
	protected:
		v2m::range_map					*m_compressed_ranges{};
		std::vector <range_el_iterator>	m_iterators;
		std::set <std::size_t>			m_valid_alts;
		
	public:
		read_compressed_vh_delegate(class generate_context &ctx, v2m::range_map &compressed_ranges):
			vh_sequence_writer(*this, ctx.reference(), ctx.null_allele_seq()),
			vh_delegate(ctx),
			m_compressed_ranges(&compressed_ranges),
			m_iterators(m_compressed_ranges->size())
		{
		}
		
		virtual bool is_valid_alt(std::size_t const alt_idx) const override { return true; }
		virtual std::set <std::size_t> const &valid_alts(v2m::variant &var) const override { return m_valid_alts; }
		virtual void prepare(v2m::vcf_reader &reader) override;
		virtual void handle_variant(v2m::variant &var) override;
		void enumerate_genotype(
			v2m::variant &var,
			std::size_t const sample_no,
			std::function <void(uint8_t, std::size_t, bool)> const &cb
		) override;
			
	protected:
		bool update_iterator_position(v2m::variant const &var, range_el_iterator &it);
	};
	
	
	class all_genotypes_vh_delegate : public vh_stats <true>, public vh_sequence_writer, public vh_delegate
	{
	public:
		all_genotypes_vh_delegate(class generate_context &ctx):
			vh_sequence_writer(*this, ctx.reference(), ctx.null_allele_seq()),
			vh_delegate(ctx)
		{
		}
		
		virtual void prepare(v2m::vcf_reader &reader) override;
		virtual void handle_variant(v2m::variant &var) override;
	};
	
	
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
	
	
	v2m::haplotype_map::iterator create_haplotype(
		v2m::haplotype_map &haplotypes,
		std::size_t const sample_no,
		std::size_t const ploidy
	)
	{
		return haplotypes.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(sample_no),
			std::forward_as_tuple(ploidy)
		).first;
	}
	
	
	v2m::haplotype_map::iterator find_or_create_haplotype(
		v2m::haplotype_map &haplotypes,
		std::size_t const sample_no,
		std::size_t const ploidy
	)
	{
		// Since file_ostream is not movable, check first if the vector has already been created.
		// If not, create it with the exact size to avoid resizing and thus moving later.
		auto it(haplotypes.find(sample_no));
		if (haplotypes.end() == it)
			it = create_haplotype(haplotypes, sample_no, ploidy);
		
		return it;
	}
	
	
	void generate_context::finish_init(
		char const *out_reference_fname,
		bool const should_reduce_samples,
		std::size_t const variant_padding,
		bool const allow_switch_to_ref
	)
	{
		if (out_reference_fname)
			m_out_reference_fname.emplace(out_reference_fname);
		
		if (should_reduce_samples)
			m_genotype_delegate.reset(new compressed_genotypes_handling_delegate(out_reference_fname != nullptr, variant_padding, allow_switch_to_ref));
		else
			m_genotype_delegate.reset(new all_genotypes_handling_delegate);
	
		m_genotype_delegate->set_generate_context(*this);
		m_variant_handler.get_variant_buffer().set_delegate(m_variant_handler);
	}
	
	
	void generate_context::check_ploidy()
	{
		size_t i(0);
		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(v2m::vcf_field::ALL);
		
		m_vcf_reader.fill_buffer();
		if (!m_vcf_reader.parse([this](v2m::transient_variant const &var) -> bool {
			for (auto const &kv : m_vcf_reader.sample_names())
			{
				auto const sample_no(kv.second);
				auto const &sample(var.sample(sample_no));
				m_ploidy[sample_no] = sample.ploidy();
			}
			
			return false;
		}))
		{
			v2m::fail("Unable to read the first variant");
		}
	}
	
	
	void generate_context::check_ref()
	{
		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(v2m::vcf_field::REF);
		bool found_mismatch(false);
		std::size_t i(0);
		
		bool should_continue(false);
		do {
			m_vcf_reader.fill_buffer();
			should_continue = m_vcf_reader.parse(
				[this, &found_mismatch, &i]
				(v2m::transient_variant const &var)
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
	
	
	void generate_context::prepare_sample_names_and_generate_sequences()
	{
		// Prepare sample names for enumeration and calculate rounds.
		{
			m_genotype_delegate->prepare_sample_names();
	
			// FIXME: the delegate's m_compressed_ranges is the wrong object.
			auto const sample_count(m_genotype_delegate->sample_count());
			m_total_rounds = std::ceil(1.0 * sample_count / m_chunk_size);
		}
		
		std::cerr << "Generating haplotype sequences…" << std::endl;
		m_start_time = std::chrono::system_clock::now();
		generate_sequences(m_out_reference_fname.operator bool());
	}
	
	
	// Handle m_chunk_size samples.
	void generate_context::generate_sequences(bool const output_reference)
	{
		// Open files for the samples. If no files were opened, exit.
		m_haplotypes.clear();
		m_genotype_delegate->update_haplotypes(
			m_haplotypes,
			output_reference
		);
		
		auto const haplotype_count(m_haplotypes.size());
		if (0 == haplotype_count)
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
		m_variant_handler.process_variants();
		// Continue in m_variant_handler_delegate's finish().
	}
	
	
	void generate_context::finish_round()
	{
		// Save the stream state.
		boost::io::ios_flags_saver ifs(std::cerr);
	
		// Change FP notation.
		std::cerr << std::fixed;
		
		auto const end_time(std::chrono::system_clock::now());
		std::chrono::duration <double> elapsed_seconds(end_time - m_round_start_time);
		std::cerr << "Finished in " << (elapsed_seconds.count() / 60.0) << " minutes." << std::endl;
	}
	
	
	void generate_context::set_variant_handler_delegate(std::unique_ptr <v2m::variant_handler_delegate> &&delegate)
	{
		m_variant_handler_delegate = std::move(delegate);
		m_variant_handler.set_delegate(*m_variant_handler_delegate);
	}
	
	
	void generate_context::load_and_generate(
		char const *reference_fname,
		char const *variants_fname,
		char const *report_fname,
		bool const should_check_ref
	)
	{
		// Open the files.
		std::cerr << "Opening files…" << std::endl;
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
				m_sv_handling_method,
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
		
		m_genotype_delegate->process_first_phase();
	}
	
	
	void genotype_handling_delegate::update_haplotypes_with_ref(v2m::haplotype_map &haplotypes)
	{
		assert(m_generate_context->has_out_reference_fname());
		
		// Check if reference output was requested.
		v2m::always_assert(
			haplotypes.cend() == haplotypes.find(v2m::REF_SAMPLE_NUMBER),
			"REF_SAMPLE_NUMBER already in use"
		);
		
		auto it(create_haplotype(haplotypes, v2m::REF_SAMPLE_NUMBER, 1));
		auto &haplotype_vec(it->second);
		open_file_for_writing(
			m_generate_context->out_reference_fname().c_str(),
			haplotype_vec[0].output_stream,
			m_generate_context->should_overwrite_files()
		);
	}
	
	
	void all_genotypes_handling_delegate::process_first_phase()
	{
		std::unique_ptr <v2m::variant_handler_delegate> delegate(new all_genotypes_vh_delegate(*this->m_generate_context));
		m_generate_context->set_variant_handler_delegate(std::move(delegate));
		m_generate_context->prepare_sample_names_and_generate_sequences();
	}
	
	
	void all_genotypes_handling_delegate::prepare_sample_names()
	{
		auto const &sample_names(m_generate_context->vcf_reader().sample_names());
		m_sample_names_it = sample_names.cbegin();
		m_sample_names_end = sample_names.cend();
		m_sample_count = sample_names.size();
	}
	
	
	void all_genotypes_handling_delegate::update_haplotypes(
		v2m::haplotype_map &haplotypes,
		bool const output_reference
	)
	{
		auto const &ploidy(m_generate_context->ploidy());
		auto const chunk_size(m_generate_context->chunk_size());
		auto const should_overwrite_files(m_generate_context->should_overwrite_files());
		
		size_t i(0);
		while (m_sample_names_it != m_sample_names_end)
		{
			auto const &sample_name(m_sample_names_it->first);
			auto const sample_no(m_sample_names_it->second);
			auto const current_ploidy(ploidy.find(sample_no)->second);
		
			auto it(find_or_create_haplotype(haplotypes, sample_no, current_ploidy));
			auto &haplotype_vec(it->second);
		
			for (size_t j(1); j <= current_ploidy; ++j)
			{
				auto const fname(boost::str(boost::format("%s-%u") % sample_name % j));
				open_file_for_writing(fname.c_str(), haplotype_vec[j - 1].output_stream, should_overwrite_files);
			}
	
			++m_sample_names_it;
			++i;
		
			if (i == chunk_size - (output_reference ? 1 : 0))
				break;
		}
		
		if (output_reference)
			update_haplotypes_with_ref(haplotypes);
	}
	
	
	void compressed_genotypes_handling_delegate::process_first_phase()
	{
		std::unique_ptr <v2m::variant_handler_delegate> delegate(
			new compress_vh_delegate(
				*m_generate_context,
				m_compressed_ranges,
				m_variant_padding,
				m_generate_context->has_out_reference_fname(),
				m_allow_switch_to_ref
			)
		);
		m_generate_context->set_variant_handler_delegate(std::move(delegate));
		m_generate_context->process_variants();
		// Continue in compress_vh_delegate's finish().
	}
	
	
	void compressed_genotypes_handling_delegate::update_haplotypes(
		v2m::haplotype_map &haplotypes,
		bool const output_reference
	)
	{
		auto const chunk_size(m_generate_context->chunk_size());
		// Check whether the reference is supposed to be output on any round.
		auto const limit(m_compressed_ranges.size() - (m_output_reference ? 1 : 0));
		
		for (std::size_t i(0); i < chunk_size; ++i)
		{
			if (limit == m_sample_idx)
				break;
			
			auto const sample_id(1 + m_sample_idx);
			auto it(find_or_create_haplotype(haplotypes, sample_id, 1));
			auto &haplotype_vec(it->second);
			auto const fname(boost::str(boost::format("%u") % sample_id));
			open_file_for_writing(fname.c_str(), haplotype_vec[0].output_stream, m_generate_context->should_overwrite_files());
			
			++m_sample_idx;
		}
		
		// Check whether the reference is supposed to be output on this round.
		if (output_reference)
			update_haplotypes_with_ref(haplotypes);
	}
	
	
	void vh_stats <true>::assigned_alt_to_sequence(std::size_t const alt_idx)
	{
		++m_non_ref_totals.handled_count;
		++m_counts_by_alt[alt_idx].handled_count;
	}
	
	
	void vh_stats <true>::found_overlapping_alt(
		std::size_t const lineno,
		uint8_t const alt_idx,
		std::size_t const sample_no,
		uint8_t const chr_idx
	)
	{
		if (m_overlapping_alts.insert(lineno).second)
		{
			std::cerr << "Overlapping alternatives on line " << lineno
			<< " for sample " << sample_no << ':' << (int) chr_idx
			<< " (and possibly others); skipping when needed." << std::endl;
		}
	
		if (this->generate_context().error_logger().is_logging_errors())
			m_skipped_samples.emplace_back(sample_no, alt_idx, chr_idx);
	}
	
	
	void vh_stats <true>::handle_variant(v2m::variant &var)
	{
		m_skipped_samples.clear();
		m_counts_by_alt.clear();
		m_non_ref_totals.reset();
	}
	
	
	void vh_stats <true>::handled_alt(std::size_t const alt_idx)
	{
		++m_non_ref_totals.total_count;
		++m_counts_by_alt[alt_idx].total_count;
	}
	
	
	void vh_stats <true>::handled_haplotypes(v2m::variant &var)
	{
		// Report errors if needed.
		auto &error_logger(this->generate_context().error_logger());
		if (error_logger.is_logging_errors())
		{
			auto const lineno(var.lineno());
			for (auto const &s : m_skipped_samples)
				error_logger.log_overlapping_alternative(lineno, s.sample_no, s.chr_idx, m_counts_by_alt[s.alt_idx], m_non_ref_totals);
		}
	}
	
	
	v2m::vcf_reader::sample_name_map const &vh_delegate::sample_names() const
	{
		return m_ctx->vcf_reader().sample_names();
	}
	
	
	bool vh_delegate::is_valid_alt(std::size_t const alt_idx) const
	{
		return m_ctx->variant_handler().is_valid_alt(alt_idx);
	}
	
	
	std::set <std::size_t> const &vh_delegate::valid_alts(v2m::variant &var) const
	{
		return m_ctx->variant_handler().valid_alts();
	}
	
	
	void vh_delegate::enumerate_genotype(
		v2m::variant &var,
		std::size_t const sample_no,
		std::function <void(uint8_t, std::size_t, bool)> const &cb
	)
	{
		return m_ctx->variant_handler().enumerate_genotype(var, sample_no, cb);
	}
	
	
	void compress_vh_delegate::prepare(v2m::vcf_reader &reader)
	{
		reader.set_parsed_fields(v2m::vcf_field::ALL);
		m_sample_reducer.prepare();
	}
	
	
	void compress_vh_delegate::handle_variant(v2m::variant &var)
	{
		vh_stats <true>::handle_variant(var);
		m_sample_reducer.handle_variant(var);
	}


	void compress_vh_delegate::finish()
	{
		m_sample_reducer.finish();

		auto &ctx(*m_ctx);
		std::unique_ptr <v2m::variant_handler_delegate> delegate(new read_compressed_vh_delegate(ctx, m_sample_reducer.compressed_ranges()));
		ctx.set_variant_handler_delegate(std::move(delegate));
		// *this is now invalid b.c. set_variant_handler_delegate replaced the unique_ptr contents that held it.
		ctx.prepare_sample_names_and_generate_sequences();
	}
	
	
	void vh_sequence_writer::finish()
	{
		m_sequence_writer.finish();
		auto &ctx(this->generate_context());
		ctx.finish_round();
		ctx.generate_sequences();
	}
	
	
	void read_compressed_vh_delegate::prepare(v2m::vcf_reader &reader)
	{
		// Fill the iterator vector.
		m_iterators.resize(m_compressed_ranges->size());
		for (auto const &ref : boost::combine(*m_compressed_ranges, m_iterators))
		{
			auto const &set(ref.get <0>());
			auto &it(ref.get <1>());
			range_el_iterator new_it(set.cbegin(), set.cend());
			it = std::move(new_it);
		}
		
		reader.set_parsed_fields(v2m::vcf_field::ALT);
		m_sequence_writer.prepare(m_ctx->haplotypes());
	}
	
	
	bool read_compressed_vh_delegate::update_iterator_position(v2m::variant const &var, range_el_iterator &it)
	{
		auto const var_pos(var.pos() - 1);
		while (!it.empty())
		{
			auto &kv(it.front()); // POS (zero-based) → variant_sequence
			auto const seq_pos(kv.first);
			
			// Check if the sequence occurs after the current variant.
			if (var_pos < seq_pos)
				return false;
			
			// Sequence beginning is before the current variant position, i.e. seq_pos ≤ var_pos
			auto const seq_end_pos(kv.second.end_pos());
			if (var_pos < seq_end_pos)
				return true;
			
			// Now seq_end_pos ≤ var_pos, check the next sequence.
			it.drop_front();
		}
		
		// This could be changed to indicate that the iterator / range may be removed from m_iterators.
		return false;
	}
	
	
	void read_compressed_vh_delegate::handle_variant(v2m::variant &var)
	{
		// For each sample, find the current range.
		m_valid_alts.clear();
		for (auto &it : m_iterators)
		{
			if (update_iterator_position(var, it))
			{
				uint8_t alt_idx(0);
				if (it.front().second.get_alt(var.lineno(), alt_idx))
					m_valid_alts.emplace(alt_idx);
			}
		}
		
		m_sequence_writer.handle_variant(var);
	}
	
	
	void read_compressed_vh_delegate::enumerate_genotype(
		v2m::variant &var,
		std::size_t const sample_no,
		std::function <void(uint8_t, std::size_t, bool)> const &cb
	)
	{
		auto const &sample_map(m_compressed_ranges->at(sample_no));
		if (0 == sample_map.size())
		{
			cb(0, 0, true);
			return;
		}

		// Find the variant_sequence that starts after the current position.
		// If the found sequence is the first one, output REF.
		auto const pos(var.pos());
		auto it(sample_map.upper_bound(pos));
		if (sample_map.cbegin() == it)
		{
			cb(0, 0, true);
			return;
		}

		// Make it point to the variant_sequence that starts before the current position.
		--it;
		uint8_t alt_idx(0);
		auto const lineno(var.lineno());
		auto const &var_seq(it->second);
		var_seq.get_alt(lineno, alt_idx);	// alt_idx remains zero if get_alt returns false.
		cb(0, alt_idx, true);
	}
	
	
	void all_genotypes_vh_delegate::prepare(v2m::vcf_reader &reader)
	{
		reader.set_parsed_fields(v2m::vcf_field::ALL);
		m_sequence_writer.prepare(m_ctx->haplotypes());
	}
	
	
	void all_genotypes_vh_delegate::handle_variant(v2m::variant &var)
	{
		vh_stats <true>::handle_variant(var);
		m_sequence_writer.handle_variant(var);
	}
}


namespace vcf2multialign {
	
	void generate_haplotypes(
		char const *reference_fname,
		char const *variants_fname,
		char const *out_reference_fname,
		char const *report_fname,
		char const *null_allele_seq,
		std::size_t const chunk_size,
		std::size_t const variant_padding,
		sv_handling const sv_handling_method,
		bool const should_overwrite_files,
		bool const should_check_ref,
		bool const should_reduce_samples,
		bool const allow_switch_to_ref
	)
	{
		dispatch_ptr <dispatch_queue_t> main_queue(dispatch_get_main_queue(), true);
		dispatch_ptr <dispatch_queue_t> parsing_queue(
			dispatch_queue_create("fi.iki.tsnorri.vcf2multialign.parsing_queue", DISPATCH_QUEUE_SERIAL),
			false
		);
		
		// generate_context needs to be allocated on the heap because later dispatch_main is called.
		// The class deallocates itself in cleanup().
		generate_context *ctx(new generate_context(
			std::move(main_queue),
			std::move(parsing_queue),
			out_reference_fname,
			null_allele_seq,
			sv_handling_method,
			chunk_size,
			variant_padding,
			should_overwrite_files,
			should_reduce_samples,
			allow_switch_to_ref
		));
			
		ctx->load_and_generate(
			reference_fname,
			variants_fname,
			report_fname,
			should_check_ref
		);
	}
}
