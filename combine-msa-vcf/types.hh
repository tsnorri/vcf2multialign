/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VCF_TYPES_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VCF_TYPES_HH

#include <libbio/vcf/variant.hh>
#include <stdexcept>
#include <string>
#include <utility>


namespace vcf2multialign {
	
	enum class segment_type : std::uint8_t
	{
		MATCH = 0,
		MISMATCH,
		DELETION,
		INSERTION,
		INSERTION_WITH_SNP,
		MIXED,
		MIXED_ALT_STARTS_WITH_GAP
	};
	
	
	enum class msa_variant_output : std::uint8_t
	{
		NONE = 0,
		ALL,
		ALT_MATCHES_REF
	};
	
	
	struct variant_record
	{
		libbio::vcf::variant	variant;
		std::size_t				aligned_position{};				// 1-based, 0 indicates invalid record.
		std::size_t				size{};
		bool					is_skipped{};
	};
	
	typedef std::vector <variant_record> variant_record_vector;
	
	
	struct sequence_character
	{
		std::size_t	position{};		// Position within the sequence.
		char		character{};
		
		sequence_character() = default;
		
		sequence_character(char const character_, std::size_t const position_):
			position(position_),
			character(std::toupper(character_))
		{
		}
	};
	
	
	struct aligned_character_pack
	{
		sequence_character	ref;
		sequence_character	alt;
		std::size_t			aligned_position{};
		char				next_ref{};
		char				next_alt{};
		
		aligned_character_pack(std::size_t const aligned_position_, sequence_character const &ref_, sequence_character const &alt_, char const next_ref_, char const next_alt_):
			ref(ref_),
			alt(alt_),
			aligned_position(aligned_position_),
			next_ref(next_ref_),
			next_alt(next_alt_)
		{
		}
	};
	
	
	struct aligned_string
	{
		std::string			string;
		std::size_t			position{};
		
		aligned_string() = default;
		
		aligned_string(std::string string_, std::size_t const position_):
			string(std::move(string_)),
			position(position_)
		{
		}
	};
	
	
	struct aligned_segment
	{
		aligned_string		ref;
		aligned_string		alt;
		std::size_t			aligned_position{};
		segment_type		type{};
		
		aligned_segment() = default;
		
		aligned_segment(
			std::string ref_,
			std::string alt_,
			std::size_t const ref_pos,
			std::size_t const alt_pos,
			std::size_t const aln_pos,
			segment_type const type_
		):
			ref(std::move(ref_), ref_pos),
			alt(std::move(alt_), alt_pos),
			aligned_position(aln_pos),
			type(type_)
		{
		}
		
		inline void reset(aligned_character_pack const &pack, segment_type const st);
		inline std::size_t alt_end() const;
	};
	
	typedef std::vector <aligned_segment> aligned_segment_vector;
	
	
	// Wrap characters into sequence_characters, count them to determine positions.
	struct position_counter
	{
		std::size_t position{};
		
		inline sequence_character add_chr(char const c);
	};
	
	
	sequence_character position_counter::add_chr(char const c)
	{
		if ('-' != c)
			++position;
		
		return sequence_character(c, position - 1);
	}
	
	
	void aligned_segment::reset(aligned_character_pack const &pack, segment_type const st)
	{
		ref.string.clear();
		alt.string.clear();
		ref.position = pack.ref.position;
		alt.position = pack.alt.position;
		aligned_position = pack.aligned_position;
		type = st;
	}
		
		
	std::size_t aligned_segment::alt_end() const
	{
		switch (type)
		{
			case segment_type::MATCH:
				return alt.position + ref.string.size();

			case segment_type::DELETION:
				return alt.position + 1;

			case segment_type::MISMATCH:
			case segment_type::INSERTION:
			case segment_type::INSERTION_WITH_SNP:
			case segment_type::MIXED:
			case segment_type::MIXED_ALT_STARTS_WITH_GAP:
				return alt.position + alt.string.size();
		}
	}
	
	
	inline std::ostream &operator<<(std::ostream &os, segment_type const st)
	{
		switch (st)
		{
			case segment_type::MATCH:
				os << "MATCH";
				return os;
				
			case segment_type::MISMATCH:
				os << "MISMATCH";
				return os;

			case segment_type::DELETION:
				os << "DELETION";
				return os;

			case segment_type::INSERTION:
				os << "INSERTION";
				return os;

			case segment_type::INSERTION_WITH_SNP:
				os << "INSERTION_WITH_SNP";
				return os;

			case segment_type::MIXED:
				os << "MIXED";
				return os;
				
			case segment_type::MIXED_ALT_STARTS_WITH_GAP:
				os << "MIXED_ALT_STARTS_WITH_GAP";
				return os;
		};
	}
	
		
	inline std::ostream &operator<<(std::ostream &os, sequence_character const &sc)
	{
		os << '[' << sc.position << "] '" << sc.character << '\'';
		return os;
	}
	
	
	inline std::ostream &operator<<(std::ostream &os, aligned_character_pack const &pack)
	{
		os << '[' << pack.aligned_position << "] Ref: " << pack.ref << " Alt: " << pack.alt << " Next_ref: '" << pack.next_ref << "' Next_alt: '" << pack.next_alt << '\'';
		return os;
	}
	
	
	inline std::ostream &operator<<(std::ostream &os, aligned_string const &str)
	{
		os << '[' << str.position << "] '" << str.string << '\'';
		return os;
	}

	
	inline std::ostream &operator<<(std::ostream &os, aligned_segment const &seg)
	{
		os << '[' << seg.aligned_position << "] Ref: " << seg.ref << " Alt: " << seg.alt << " Type: " << seg.type;
		return os;
	}
}

#endif
