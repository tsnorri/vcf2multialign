/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch_fn.hh>
#include <unistd.h>
#include <vcf2multialign/generate_haplotypes.hh>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	v2m::sv_handling sv_handling_method(enum_structural_variants const sva)
	{
		switch (sva)
		{
			case structural_variants_arg_discard:
				return v2m::sv_handling::DISCARD;
				
			case structural_variants_arg_keep:
				return v2m::sv_handling::KEEP;
				
			case structural_variants__NULL:
			default:
				lb::fail("Unexpected value for structural variant handling.");
				return v2m::sv_handling::KEEP; // Not reached.
		}
	}
}


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

	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.

	// libdispatch on macOS does not need pthread_workqueue.
#ifdef __linux__
	pthread_workqueue_init_np();
#endif

	v2m::generate_haplotypes(
		args_info.reference_arg,
		args_info.variants_arg,
		args_info.output_reference_arg,
		args_info.report_file_arg,
		args_info.null_allele_seq_arg,
		args_info.chunk_size_arg,
		sv_handling_method(args_info.structural_variants_arg),
		args_info.overwrite_flag,
		!args_info.no_check_ref_flag
	);
		
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	return EXIT_SUCCESS;
}
