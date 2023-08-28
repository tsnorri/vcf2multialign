/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TRANSPOSE_MATRIX_HH
#define VCF2MULTIALIGN_TRANSPOSE_MATRIX_HH

#include <libbio/int_matrix.hh>


namespace vcf2multialign {
	
	libbio::bit_matrix transpose_matrix(libbio::bit_matrix const &mat);
}

#endif
