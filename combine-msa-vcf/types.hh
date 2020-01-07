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
		MIXED
	};
	
	
	// From C++20.
	template <typename t_type>
	struct remove_cvref
	{
		typedef std::remove_cv_t <std::remove_reference_t <t_type>> type;
	};
	
	template <typename t_type>
	using remove_cvref_t = typename remove_cvref <t_type>::type;
	
	
	struct variant_record
	{
		libbio::variant	variant;
		std::size_t		aligned_position{};				// 1-based, 0 indicates invalid record.
		std::size_t		size{};
		bool			is_accounted_for_overlaps{};	// Bookkeeping.
	};
	
	
	struct sequence_character
	{
		std::size_t	position{};		// Position within the sequence.
		char		character{};
		
		sequence_character() = default;
		
		sequence_character(char const character_, std::size_t const position_):
			position(position_),
			character(character_)
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
	};
	
	
	struct aligned_segment
	{
		aligned_string		ref;
		aligned_string		alt;
		std::size_t			aligned_position{};
		segment_type		type{};
		
		inline void reset(aligned_character_pack const &pack, segment_type const st);
	};
	
	
	// Wrap characters into sequence_characters, count them to determine positions.
	struct position_counter
	{
		std::size_t position{};
		
		inline sequence_character add_chr(char const c);
	};
	
	
	// May be used in place of an output iterator to forward arguments of operator=.
	template <typename t_dst>
	class forwarder
	{
	public:
		typedef std::ptrdiff_t difference_type; // Needed for range-v3 only.
		
	protected:
		t_dst	*m_dst{};
	
	public:
		// ranges::semiregular<T> (in range-v3/include/concepts/concepts.hpp) requires copyable and default_constructible.
		// m_dst is required, though, so throw in case the default constructor is somehow called.
		forwarder()
		{
			throw std::runtime_error("forwarder’s default constructor should not be called.");
		}
		
		forwarder(t_dst &dst):
			m_dst(&dst)
		{
		}
	
		forwarder &operator++() { return *this; }				// Return *this.
		forwarder &operator++(int) { return *this; }			// Return *this.
		forwarder &operator*() { return *this; }				// Return *this.
		forwarder const &operator*() const { return *this; }	// Return *this.
	
		template <typename t_arg>
		forwarder &operator=(t_arg &&arg) { m_dst->handle(std::forward <t_arg>(arg)); return *this; }
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
