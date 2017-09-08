/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TYPES_HH
#define VCF2MULTIALIGN_TYPES_HH

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <set>
#include <vector>


namespace vcf2multialign {
	typedef std::vector <char> vector_type;
	typedef std::set <std::size_t> variant_set;
	
	
	struct sample_count
	{
		std::size_t handled_count{0};
		std::size_t total_count{0};
		
		void reset() { handled_count = 0; total_count = 0; }
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
	
	
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_source>	file_istream;
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_sink>	file_ostream;
	
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
