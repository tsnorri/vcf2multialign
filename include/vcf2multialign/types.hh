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
	
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_source>	file_istream;
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_sink>	file_ostream;
}

#endif
