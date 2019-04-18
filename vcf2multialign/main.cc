/*
 * Copyright (c) 2017-2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <unistd.h>
#include <vcf2multialign/generate_haplotypes.hh>
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
	
	if (args_info.chunk_size_arg <= 0)
	{
		std::cerr << "Chunk size must be positive." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	v2m::output_type ot(v2m::output_type::SEQUENCE_FILES);
	if (args_info.output_sequences_given)
		;
	else if (args_info.output_variant_graph_given)
		ot = v2m::output_type::VARIANT_GRAPH;
	else
		libbio_fail("Unexpected output type.");

	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

	v2m::generate_haplotypes(
		args_info.reference_arg,
		args_info.variants_arg,
		args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
		args_info.chromosome_given ? args_info.chromosome_arg : nullptr,
		ot,
		args_info.output_arg,
		args_info.output_reference_arg,
		args_info.report_file_arg,
		args_info.null_allele_seq_arg,
		args_info.chunk_size_arg,
		args_info.overwrite_flag,
		!args_info.no_check_ref_flag
	);
		
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
