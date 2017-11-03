/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TYPES_HH
#define VCF2MULTIALIGN_TYPES_HH

#include <map>
#include <set>
#include <vcf2multialign/file_handling.hh>
#include <vector>


namespace vcf2multialign {
	enum { REF_SAMPLE_NUMBER = 0 };
	
	typedef std::vector <char> vector_type;
	typedef std::set <std::size_t> variant_set;
	typedef std::map <char const *, std::size_t> subgraph_map;
	
	struct sample_count
	{
		std::size_t handled_count{0};
		std::size_t total_count{0};
		
		void reset() { handled_count = 0; total_count = 0; }
	};
	
	template <typename t_ostream>
	struct haplotype
	{
		size_t current_pos{0};
		t_ostream output_stream;
	};
	
	template <typename t_ostream>
	using haplotype_map = std::map <
		std::size_t,							// Sample (line) number
		std::vector <haplotype <t_ostream>>		// All haplotype sequences
	>;

	template <typename t_ostream>
	using haplotype_ptr_map = std::map <
		std::size_t,							// Sample (line) number
		std::vector <haplotype <t_ostream> *>	// Haplotype sequences by chromosome index
	>;
	
	template <typename t_ostream>
	using alt_map = std::map <
		std::string,							// ALT
		haplotype_ptr_map <t_ostream>
	>;
	
	typedef std::map <std::size_t, std::size_t> ploidy_map;
	
	enum class vcf_field : uint8_t {
		CHROM	= 0,
		POS		= 1,
		ID		= 2,
		REF		= 3,
		ALT		= 4,
		QUAL	= 5,
		FILTER	= 6,
		INFO	= 7,
		FORMAT	= 8,
		ALL		= 9,
		VARY	= 10
	};
	
	enum class format_field : uint8_t {
		GT		= 0,
		DP,
		GQ,
		PS,
		PQ,
		MQ
	};
	
	enum { NULL_ALLELE = std::numeric_limits <uint8_t>::max() };
	
	enum class sv_handling : uint8_t {
		DISCARD		= 0,
		KEEP
	};
	
	enum class sv_type : uint8_t {
		NONE		= 0,
		DEL,
		INS,
		DUP,
		INV,
		CNV,
		DUP_TANDEM,
		DEL_ME,
		INS_ME,
		UNKNOWN
	};
	
	char const *to_string(sv_type const svt);
	
	
	struct graph_range
	{
		std::size_t	range_start_offset{0};
		std::size_t	range_length{0};
		std::size_t start_lineno{0};
		std::size_t	variant_count{0};
		
		graph_range() = default;
		graph_range(
			std::size_t	range_start_offset_,
			std::size_t	range_length_,
			std::size_t start_lineno_,
			std::size_t	variant_count_
		):
			range_start_offset(range_start_offset_),
			range_length(range_length_),
			start_lineno(start_lineno_),
			variant_count(variant_count_)
		{
		}
	};
}

#endif
