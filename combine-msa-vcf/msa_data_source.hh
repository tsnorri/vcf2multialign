/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_MSA_DATA_SOURCE_HH
#define VCF2MULTIALIGN_COMBINE_MSA_MSA_DATA_SOURCE_HH

#include <vcf2multialign/types.hh>
#include "overlap_counter.hh"
#include "types.hh"
#include "vcf_record_generator.hh"


namespace vcf2multialign {
	
	struct msa_segmentation
	{
		aligned_segment_vector				overlapping_segments;
		variant_record_vector				overlapping_variants;
		class overlap_counter				overlap_counter;
	};
	
	
	struct segmentation_handler
	{
		virtual void process_variants(msa_segmentation &seg) = 0;
	};
	
	
	class msa_data_source
	{
		friend class forwarder <msa_data_source>;
		
	protected:
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
		
	protected:
		fsm									m_fsm;
		msa_segmentation					m_segmentation;
		aligned_segment						m_current_segment;
		segmentation_handler				*m_segmentation_handler{};
		std::size_t							m_handled_variants{};
		std::size_t							m_handled_characters{};
		std::uint16_t						m_ploidy{1};
		bool								m_logs_status{};
		
	public:
		msa_data_source() = default;
		
		msa_data_source(std::uint16_t const ploidy, bool const logs_status):
			m_ploidy(ploidy),
			m_logs_status(logs_status)
		{
		}
		
		void process_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen, segmentation_handler &handler);
		
	protected:
		// Handle aligned_character_packs.
		void handle(aligned_character_pack const &&pack); // The pack only contains primitive types and it’s not going to be used after calling this.
		
		// Handle variant_records.
		void handle(variant_record &&rec);
		
		gap_start_position check_gaps_at_start(vector_type const &ref, vector_type const &alt) const;
		
		void prepare_msa_parser();
		void parse_msa(aligned_character_pack const &pack);
		void push_current_segment();
		void finish_msa() { push_current_segment(); }
	};
	
	
	auto msa_data_source::fsm::operator=(fsm const &other) & -> msa_data_source::fsm &
	{
		if (this != &other)
		{
			characters = other.characters;
			p = characters.data();
			pe = p + characters.size();
		}
		return *this;
	}
	
	
	auto msa_data_source::fsm::operator=(fsm &&other) & -> msa_data_source::fsm &
	{
		return this->operator=(other);
	}
	
	
	void msa_data_source::fsm::update_characters(aligned_character_pack const &pack)
	{
		characters[0] = pack.ref.character;
		characters[1] = pack.alt.character;
		characters[2] = pack.next_ref;
		characters[3] = pack.next_alt;
		p = characters.data();
	}
}

#endif
