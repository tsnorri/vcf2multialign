/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_HANDLER_HH
#define VCF2MULTIALIGN_VARIANT_HANDLER_HH

#include <dispatch/dispatch.h>
#include <map>
#include <stack>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/vcf_reader.hh>
#include <vcf2multialign/vector_source.hh>


namespace vcf2multialign {

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
	
	
	struct variant_overlap
	{
		size_t	start_pos{0};
		size_t	current_pos{0};
		size_t	end_pos{0};
		size_t	heaviest_path_length{0};
		alt_map	alt_haplotypes;
		
		variant_overlap(size_t const s, size_t const c, size_t const e, size_t const h):
			start_pos(s),
			current_pos(c),
			end_pos(e),
			heaviest_path_length(h)
		{
			if (! (start_pos <= end_pos))
				throw std::runtime_error("Bad offset order");
		}
		
		variant_overlap(size_t const s, size_t const c, size_t const e, size_t const h, alt_map &alts):
			start_pos(s),
			current_pos(c),
			end_pos(e),
			heaviest_path_length(h),
			alt_haplotypes(std::move(alts))
		{
			if (! (start_pos <= end_pos))
				throw std::runtime_error("Bad offset order");
		}
	};

	
	class variant_handler
	{
	protected:
		typedef std::stack <variant_overlap>	overlap_stack_type;
		typedef std::vector <size_t>			sample_number_vector;

	protected:
		dispatch_queue_t								m_main_queue;
		dispatch_queue_t								m_parsing_queue;
		std::function <void(void)>						m_finish_callback;
		
		vector_type	const								*m_reference{};
		
		vcf_reader										*m_vcf_reader{};
		variant											m_var;
		overlap_stack_type								m_overlap_stack;
		
		variant_set const								*m_skipped_variants{};
		variant_set										m_overlapping_alts{};
		haplotype_ptr_map								m_ref_haplotype_ptrs;	// Haplotypes to which the reference sequence is to be output.
		std::set <size_t>								m_valid_alts;
		
		vector_source <variant::sample_field_vector>	m_sample_vs;
		vector_source <variant::genotype_vector>		m_gt_vs;
		
		haplotype_map									*m_all_haplotypes{};
		alt_map											m_alt_haplotypes;
		
		std::string const								*m_null_allele_seq{};
		std::size_t										m_i{0};
		
	public:
		variant_handler(
			dispatch_queue_t main_queue,
			dispatch_queue_t parsing_queue,
			vcf_reader &vcf_reader_,
			vector_type const &reference,
			variant_set const &skipped_variants,
			std::string const &null_allele,
			std::function <void(void)> finish_callback
		):
			m_main_queue(main_queue),
			m_parsing_queue(parsing_queue),
			m_finish_callback(finish_callback),
			m_reference(&reference),
			m_vcf_reader(&vcf_reader_),
			m_var(vcf_reader_.sample_count()),
			m_skipped_variants(&skipped_variants),
			m_null_allele_seq(&null_allele)
		{
			dispatch_retain(m_main_queue);
			dispatch_retain(m_parsing_queue);
			m_var.add_format_field("GT");
		}
		
		~variant_handler()
		{
			dispatch_release(m_main_queue);
			dispatch_release(m_parsing_queue);
		}
		
		variant_handler() = default;
		variant_handler(variant_handler const &) = delete;
		variant_handler(variant_handler &&) = default;
		variant_handler &operator=(variant_handler const &) & = delete;
		variant_handler &operator=(variant_handler &&) & = default;
		
	public:
		void process_variants(haplotype_map &haplotypes);

	protected:
		void process_next_variant();
		void fill_streams(haplotype_ptr_map &haplotypes, size_t const fill_amt) const;
		void output_reference(std::size_t const output_start_pos, std::size_t const output_end_pos);
		void process_overlap_stack(size_t const var_pos);
	};
}

#endif
