/*
 * Copyright (c) 2017-2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_CAN_HANDLE_VARIANT_ALTS_HH
#define VCF2MULTIALIGN_UTILITY_CAN_HANDLE_VARIANT_ALTS_HH

#include <libbio/vcf/variant.hh>


namespace vcf2multialign {
		
	bool can_handle_variant_alts(libbio::vcf::transient_variant const &var);
	bool can_handle_variant_alt(libbio::vcf::variant_alt_base const &alt);
}

#endif
