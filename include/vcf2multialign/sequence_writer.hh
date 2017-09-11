/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_SEQUENCE_WRITER_HH
#define VCF2MULTIALIGN_SEQUENCE_WRITER_HH

#include <map>
#include <stack>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/util.hh>
#include <vcf2multialign/variant_processor_delegate.hh>


namespace vcf2multialign {

	struct sequence_writer_delegate : public virtual variant_processor_delegate
	{
		virtual std::set <std::size_t> const &valid_alts(variant &var) const = 0;
	};
	
	
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
	
	
	class sequence_writer
	{
	protected:
		typedef std::stack <variant_overlap>			overlap_stack_type;
		typedef std::vector <size_t>					sample_number_vector;
		
	protected:
		sequence_writer_delegate						*m_delegate{};
		
		vector_type	const								*m_reference{};

		overlap_stack_type								m_overlap_stack;
		
		haplotype_ptr_map								m_ref_haplotype_ptrs;	// Haplotypes to which the reference sequence is to be output.
		
		haplotype_map									*m_all_haplotypes{};
		alt_map											m_alt_haplotypes;

		std::string const								*m_null_allele_seq{};
		
	public:
		sequence_writer(
			vector_type const &reference,
			std::string const &null_allele
		):
			m_reference(&reference),
			m_null_allele_seq(&null_allele)
		{
		}
		
		void set_delegate(sequence_writer_delegate &delegate) { m_delegate = &delegate; }
		
		void prepare(haplotype_map &haplotypes);
		void handle_variant(variant &var);
		void finish();
		
	protected:
		void fill_streams(haplotype_ptr_map &haplotypes, size_t const fill_amt) const;
		void output_reference(std::size_t const output_start_pos, std::size_t const output_end_pos);
		std::size_t process_overlap_stack(size_t const var_pos);
	};
}

#endif
