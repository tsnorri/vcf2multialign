/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_OUTPUT_HANDLER_HH
#define VCF2MULTIALIGN_OUTPUT_HANDLER_HH

#include "variant_description.hh"

namespace vcf2multialign {
	
	struct output_handler
	{
		virtual ~output_handler() {}
		virtual void handle_variant_description(variant_description &&desc) = 0;
	};
}

#endif
