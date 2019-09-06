/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_TYPES_HH
#define VCF2MULTIALIGN_PREPROCESS_TYPES_HH

#include <cstdint>
#include <vector>


namespace vcf2multialign {

	typedef std::size_t							position_type;		// FIXME: replace instances of std::size_t with this.
	typedef std::vector <position_type>			position_vector;	// FIXME: replace instances of std::vector <std::size_t> with this.
	typedef std::uint16_t						path_number_type;
	typedef std::uint32_t						sample_number_type;
	typedef std::vector <path_number_type>		path_vector;
	typedef std::vector <path_vector>			path_matrix;
	typedef std::vector <sample_number_type>	sample_vector;
	typedef std::vector <sample_vector>			sample_map;
	
	enum {
		PATH_NUMBER_MAX		= std::numeric_limits <path_number_type>::max(),
		SAMPLE_NUMBER_MAX	= std::numeric_limits <sample_number_type>::max()
	};
}

#endif
