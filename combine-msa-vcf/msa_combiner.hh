/*
 * Copyright (c) 2019–2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_MSA_COMBINER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_MSA_COMBINER_HH

#include <array>
#include <vcf2multialign/types.hh>
#include <vector>
#include "overlap_counter.hh"
#include "vcf_record_generator.hh"
#include "types.hh"


namespace vcf2multialign {
	
	class msa_combiner
	{
	public:
		typedef std::vector <aligned_segment> aligned_segment_vector;
		
	protected:
		enum class variant_origin : std::uint8_t
		{
			MSA,
			VC
		};
		
		struct gap_start_position
		{
			std::size_t	aligned_start{};				// 0-based index of the first non-gap character in either sequence.
			std::size_t	ref_pad_diff{};					// Given rp (number of gap characters in the beginning of ref) and ap (same for alt), ref_pad_diff = max(0, rp - ap).
			char		first_ref_chr_after_pad{};
			bool		ref_chr_needs_left_align{};
			
			gap_start_position(
				std::size_t	aligned_start_,
				std::size_t	ref_pad_diff_,
				char		first_ref_chr_after_pad_,
				bool		ref_chr_needs_left_align_
			):
				aligned_start(aligned_start_),
				ref_pad_diff(ref_pad_diff_),
				first_ref_chr_after_pad(first_ref_chr_after_pad_),
				ref_chr_needs_left_align(ref_chr_needs_left_align_)
			{
			}
		};
		
		struct fsm
		{
			std::array <char, 4> characters{};
			char *p{};
			char *pe{};
			int cs{};
			
			fsm():
				p(characters.data()),
				pe(p + characters.size())
			{
				
			}
			
			fsm(fsm &&other):
				characters(other.characters),
				p(characters.data()),
				pe(p + characters.size())
			{
			}
			
			inline fsm &operator=(fsm const &other) &;
			inline fsm &operator=(fsm &&other) &;
			inline void update_characters(aligned_character_pack const &pack);
		};
		
		struct variant_description
		{
			// The genotype is stored in a std::vector <bool> here.
			// A pair of integers could be used instead, though, since
			// at the moment we only handle unphased VCF files and thus
			// the order of the values does not matter.
			
			std::string			ref;
			std::string			alt;
			std::vector <bool>	genotype;
			std::size_t			position{};
			std::int32_t		overlap_count{};
			variant_origin		origin{};
			bool				is_skipped{};
			
			variant_description() = default;
			
			variant_description(std::size_t const ploidy, std::int32_t const gt_count, variant_origin const origin_):
				genotype(ploidy, true),
				origin(origin_)
			{
				libbio_assert_lte(gt_count, ploidy);
				// Assign the GT values.
				genotype.resize(gt_count);
				genotype.resize(ploidy, false);
			}
			
			template <typename t_string_1, typename t_string_2>
			variant_description(
				std::size_t const position_,
				t_string_1 &&ref_,
				t_string_2 &&alt_,
				std::size_t const ploidy,
				std::int32_t const overlap_count_,
				variant_origin const origin_
			):
				ref(std::forward <t_string_1>(ref_)),
				alt(std::forward <t_string_2>(alt_)),
				genotype(ploidy, false),
				position(position_),
				overlap_count(overlap_count_),
				origin(origin_)
			{
				libbio_assert_lte(overlap_count, ploidy);
				// Assign the GT values.
				genotype.resize(overlap_count);
				genotype.resize(ploidy, true);
			}
		};
		
	protected:
		fsm									m_fsm;
		aligned_segment_vector				m_overlapping_segments;
		std::vector <variant_record>		m_overlapping_variants;
		std::vector <variant_description>	m_output_variants;
		overlap_counter						m_overlap_counter;
		aligned_segment						m_current_segment;
		std::string							m_output_chr_id;
		std::ostream						*m_os{};
		std::size_t							m_max_rec_end{}; // Max. end co-ordinate of a variant encountered relative to the ad-hoc reference.
		std::size_t							m_handled_variants{};
		std::size_t							m_handled_characters{};
		std::uint16_t						m_ploidy{1};
		bool								m_need_alt_pos{true};
		bool								m_logs_status{};
	
	public:
		msa_combiner() = default;
		
		msa_combiner(std::string output_chr_id, std::uint16_t const ploidy, std::ostream &os, bool const logs_status):
			m_output_chr_id(std::move(output_chr_id)),
			m_os(&os),
			m_ploidy(ploidy),
			m_logs_status(logs_status)
		{
		}
		
		// Handle aligned_character_packs.
		void handle(aligned_character_pack const &&pack); // The pack only contains primitive types and it’s not going to be used after calling this.
		
		// Handle variant_records.
		void handle(variant_record &&rec);
		
		// Entry point.
		void process_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen);
		
	protected:
		void push_current_segment();
		void merge_output_variants(std::size_t const partition_point);
		
		void handle_one_segment_msa(
			aligned_segment const &seg,
			std::int32_t overlap_count,
			overlap_counter::const_iterator overlap_it,
			overlap_counter::const_iterator const overlap_end
		);
		
		std::size_t process_variants_msa(
			std::size_t const max_alt_pos,
			aligned_segment_vector::const_iterator seg_it,
			aligned_segment_vector::const_iterator const seg_end
		);
		
		void process_variant_in_range(
			variant_record const &var,
			aligned_segment_vector::const_iterator aligned_segment_begin,
			aligned_segment_vector::const_iterator const aligned_segment_end,
			std::int32_t const max_overlaps
		);

		void process_variants();
		void output_vcf_header() const;
		void filter_processed_variants_and_output(std::size_t const min_unhandled_ref_pos);
		gap_start_position check_gaps_at_start(vector_type const &ref, vector_type const &alt) const;
		
		void prepare_msa_parser();
		void parse_msa(aligned_character_pack const &pack);
		void finish_msa() { push_current_segment(); }
	};
	
	
	void msa_combiner::fsm::update_characters(aligned_character_pack const &pack)
	{
		characters[0] = pack.ref.character;
		characters[1] = pack.alt.character;
		characters[2] = pack.next_ref;
		characters[3] = pack.next_alt;
		p = characters.data();
	}
	
	
	auto msa_combiner::fsm::operator=(fsm const &other) & -> msa_combiner::fsm &
	{
		if (this != &other)
		{
			characters = other.characters;
			p = characters.data();
			pe = p + characters.size();
		}
		return *this;
	}
	
	
	auto msa_combiner::fsm::operator=(fsm &&other) & -> msa_combiner::fsm &
	{
		return this->operator=(other);
	}
	
	
	// Find the segments that overlap with the given variant.
	// The range from seg_it to seg_end needs to consist of consequtive segments sorted by their position.
	template <typename t_iterator, typename t_variant>
	inline auto find_overlapping_segment_range(
		t_iterator const seg_it,
		t_iterator const seg_end,
		t_variant const &var
	) -> std::pair <t_iterator, t_iterator>
	{
		struct {
			bool operator()(t_variant const &var, aligned_segment const &seg) const
			{
				// lt. if the variant is located before this segment.
				auto const var_pos(var.variant.zero_based_pos());
				auto const var_end_pos(var_pos + var.size);
				return var_end_pos <= seg.alt.position;
			}
			
			bool operator()(aligned_segment const &seg, t_variant const &var) const
			{
				// lt. if the variant start is located after this segment.
				auto const var_pos(var.variant.zero_based_pos());
				auto const seg_end_pos(seg.alt_end());
				return seg_end_pos <= var_pos;
			}
		} seg_cmp;
		
		return std::equal_range(seg_it, seg_end, var, seg_cmp);
	}
}

#endif
