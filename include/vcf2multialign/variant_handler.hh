/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_HANDLER_HH
#define VCF2MULTIALIGN_VARIANT_HANDLER_HH

#include <dispatch/dispatch.h>
#include <map>
#include <stack>
#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_buffer.hh>
#include <vcf2multialign/vcf_reader.hh>
#include <vcf2multialign/vector_source.hh>


namespace vcf2multialign {

	enum { REF_SAMPLE_NUMBER = 0 };
	
	struct haplotype
	{
		size_t current_pos{0};
		file_ostream output_stream;
	};
	
	typedef std::map <
		std::size_t,				// Sample (line) number
		std::vector <haplotype>		// All haplotype sequences
	> haplotype_map;

	typedef std::map <
		std::size_t,				// Sample (line) number
		std::vector <haplotype *>	// Haplotype sequences by chromosome index
	> haplotype_ptr_map;
	
	typedef std::map <
		std::string,				// ALT
		haplotype_ptr_map
	> alt_map;
	
	
	struct skipped_sample
	{
		std::size_t	sample_no{0};
		uint8_t		alt_idx{0};
		uint8_t		chr_idx{0};
		
		skipped_sample(std::size_t const sample_no_, uint8_t const alt_idx_, uint8_t const chr_idx_):
			sample_no(sample_no_),
			alt_idx(alt_idx_),
			chr_idx(chr_idx_)
		{
		}
	};
	
	
	struct variant_overlap
	{
		size_t	start_pos{0};
		size_t	current_pos{0};
		size_t	end_pos{0};
		size_t	heaviest_path_length{0};
		size_t	lineno{0};
		alt_map	alt_haplotypes;
		
		variant_overlap(
			size_t const start_pos_,
			size_t const current_pos_,
			size_t const end_pos_,
			size_t const heaviest_path_length_,
			size_t const lineno_
		):
			start_pos(start_pos_),
			current_pos(current_pos_),
			end_pos(end_pos_),
			heaviest_path_length(heaviest_path_length_),
			lineno(lineno_)
		{
			always_assert(start_pos <= end_pos, "Bad offset order");
		}
		
		variant_overlap(
			size_t const start_pos_,
			size_t const current_pos_,
			size_t const end_pos_,
			size_t const heaviest_path_length_,
			size_t const lineno_,
			alt_map &alts
		):
			start_pos(start_pos_),
			current_pos(current_pos_),
			end_pos(end_pos_),
			heaviest_path_length(heaviest_path_length_),
			lineno(lineno_),
			alt_haplotypes(std::move(alts))
		{
			always_assert(start_pos <= end_pos, "Bad offset order");
		}
	};

	
	class variant_handler : public variant_buffer_delegate
	{
	protected:
		typedef std::stack <variant_overlap>			overlap_stack_type;
		typedef std::vector <size_t>					sample_number_vector;

	protected:
		dispatch_ptr <dispatch_queue_t>					m_main_queue{};
		dispatch_ptr <dispatch_queue_t>					m_parsing_queue{};
		std::function <void(void)>						m_finish_callback;
		
		error_logger									*m_error_logger{};
		
		vector_type	const								*m_reference{};
		
		variant_buffer									m_variant_buffer;
		overlap_stack_type								m_overlap_stack;
		
		variant_set const								*m_skipped_variants{};
		variant_set										m_overlapping_alts{};
		haplotype_ptr_map								m_ref_haplotype_ptrs;			// Haplotypes to which the reference sequence is to be output.
		std::set <size_t>								m_valid_alts;
		std::vector <skipped_sample>					m_skipped_samples;				// In current variant.
		std::map <uint8_t, sample_count>				m_counts_by_alt;				// In current variant.
		sample_count									m_non_ref_totals;				// In current variant.
		
		haplotype_map									*m_all_haplotypes{};
		alt_map											m_alt_haplotypes;
		
		std::string const								*m_null_allele_seq{};
		sv_handling										m_sv_handling_method{};
		std::size_t										m_i{0};
		
	public:
		variant_handler(
			dispatch_ptr <dispatch_queue_t> const &main_queue,
			dispatch_ptr <dispatch_queue_t> const &parsing_queue,
			vcf_reader &vcf_reader_,
			vector_type const &reference,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants,
			std::string const &null_allele,
			error_logger &error_logger,
			std::function <void(void)> finish_callback
		):
			m_main_queue(main_queue),
			m_parsing_queue(parsing_queue),
			m_finish_callback(finish_callback),
			m_error_logger(&error_logger),
			m_reference(&reference),
			m_variant_buffer(vcf_reader_, main_queue, *this),
			m_skipped_variants(&skipped_variants),
			m_null_allele_seq(&null_allele),
			m_sv_handling_method(sv_handling_method)
		{
		}
		
		variant_handler() = default;
		
	public:
		void process_variants(haplotype_map &haplotypes);
		variant_buffer &get_variant_buffer() { return m_variant_buffer; }

	protected:
		virtual void handle_variant(variant &var);
		virtual void finish();

		void fill_streams(haplotype_ptr_map &haplotypes, size_t const fill_amt) const;
		void output_reference(std::size_t const output_start_pos, std::size_t const output_end_pos);
		std::size_t process_overlap_stack(size_t const var_pos);
	};
}

#endif
