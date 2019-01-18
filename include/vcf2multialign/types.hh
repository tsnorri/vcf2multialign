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
	enum { REF_SAMPLE_NUMBER = 0 };
	
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_source>	file_istream;
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_sink>	file_ostream;
	
	typedef std::vector <char> vector_type;
	typedef std::set <std::size_t> variant_set;
	
	enum class output_type : uint8_t {
		SEQUENCE_FILES	= 0,
		VARIANT_GRAPH
	};
	
	enum class sv_handling : uint8_t {
		DISCARD		= 0,
		KEEP
	};
	
	struct sample_count
	{
		std::size_t handled_count{0};
		std::size_t total_count{0};
		
		void reset() { handled_count = 0; total_count = 0; }
	};
}

#endif
