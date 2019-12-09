/*
 * Copyright (c) 2017-2018 Tuukka Norri
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
	typedef std::vector <char> vector_type;
	typedef std::set <std::size_t> variant_set;
	
	enum class output_type : std::uint8_t {
		SEQUENCE_FILES	= 0,
		VARIANT_GRAPH
	};
	
	struct sample_count
	{
		std::size_t handled_count{0};
		std::size_t total_count{0};
		
		void reset() { handled_count = 0; total_count = 0; }
	};
	
	typedef std::map <std::size_t, std::size_t> ploidy_map;
}

#endif
