/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TYPES_HH
#define VCF2MULTIALIGN_TYPES_HH

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <map>
#include <set>
#include <vector>


namespace vcf2multialign {
	enum { REF_SAMPLE_NUMBER = 0 };
	
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_source>	file_istream;
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_sink>	file_ostream;
	
	typedef std::vector <char> vector_type;
	typedef std::set <std::size_t> variant_set;
	
	struct sample_count
	{
		std::size_t handled_count{0};
		std::size_t total_count{0};
		
		void reset() { handled_count = 0; total_count = 0; }
	};
	
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
}

#endif
