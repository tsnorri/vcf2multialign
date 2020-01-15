/*
 * Copyright (c) 2019–2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_MSA_COMBINER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_MSA_COMBINER_HH

#include <array>
#include <vcf2multialign/types.hh>
#include <vector>
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
		
		struct overlap_count
		{
			std::size_t		position{};
			std::int32_t	running_sum{};
			std::int8_t		count{};
			
			overlap_count() = default;
			overlap_count(std::size_t const position_, std::int8_t const count_):
				position(position_),
				count(count_)
			{
			}
		};
		
		friend std::ostream &operator<<(std::ostream &, overlap_count const &);
		typedef std::vector <overlap_count> overlap_count_vector;
		
	protected:
		fsm									m_fsm;
		aligned_segment_vector				m_overlapping_segments;
		std::vector <variant_record>		m_overlapping_variants;
		std::vector <variant_description>	m_output_variants;
		overlap_count_vector				m_overlap_counts;
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
		void handle_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen);
		
	protected:
		void output_vcf_header() const;
		gap_start_position check_gaps_at_start(vector_type const &ref, vector_type const &alt) const;
		std::int32_t count_set_genotype_values(libbio::variant const &var, std::uint16_t const alt_idx) const;
		void push_current_segment();
		inline aligned_segment const &find_segment_for_alt_position(std::size_t const pos) const;
		void handle_overlaps(std::size_t const max_alt_pos);
		std::vector <variant_record>::iterator partition_overlapping_variants_by_pos(std::size_t const max_alt_pos);
		std::size_t handle_overlaps_msa(
			std::size_t const max_alt_pos,
			aligned_segment_vector::const_iterator begin,
			aligned_segment_vector::const_iterator const end
		);
		void handle_one_segment_msa(
			aligned_segment const &seg,
			std::int32_t const overlap_count,
			overlap_count_vector::iterator overlap_it,
			overlap_count_vector::iterator const overlap_end
		);

		void process_variant_in_range(
			variant_record const &var,
			aligned_segment_vector::const_iterator begin,
			aligned_segment_vector::const_iterator const end,
			std::int32_t const max_overlaps
		);
		void filter_processed_variants_and_output();
		
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
	
	
	inline std::ostream &operator<<(std::ostream &os, msa_combiner::overlap_count const &oc)
	{
		os << "pos: " << oc.position << " rs: " << +oc.running_sum << " count: " << +oc.count;
		return os;
	}
	
	
	auto msa_combiner::find_segment_for_alt_position(std::size_t const pos) const -> aligned_segment const &
	{
		libbio_always_assert(m_current_segment.alt.position <= pos);
		return m_current_segment;
	}
}

#endif
