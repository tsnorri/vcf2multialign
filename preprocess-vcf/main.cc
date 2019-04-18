/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <unistd.h>
#include <vcf2multialign/preprocess/preprocess_variants.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


int main(int argc, char **argv)
{
#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif

	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

	std::vector <std::string> field_names_for_filter_if_set(args_info.filter_fields_set_given);
	for (std::size_t i(0); i < args_info.filter_fields_set_given; ++i)
		field_names_for_filter_if_set[i] = args_info.filter_fields_set_arg[i];
	
	v2m::preprocess_variants(
		args_info.reference_arg,
		args_info.variants_arg,
		args_info.output_variants_arg,
		args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
		args_info.chromosome_given ? args_info.chromosome_arg : nullptr,
		field_names_for_filter_if_set,
		args_info.overwrite_flag
	);
	
	return EXIT_SUCCESS;
}
