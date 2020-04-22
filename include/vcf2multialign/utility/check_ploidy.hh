/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_CHECK_PLOIDY_HH
#define VCF2MULTIALIGN_UTILITY_CHECK_PLOIDY_HH

#include <cstddef>
#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {

	void check_ploidy(libbio::vcf::reader &vcf_reader, ploidy_map &out_ploidy);
}

#endif
