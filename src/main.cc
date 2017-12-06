/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/generate_haplotypes.hh>
#include <vcf2multialign/gzip_sink.hh>
#include <vcf2multialign/util.hh>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


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
				v2m::fail("Unexpected value for structural variant handling.");
				return v2m::sv_handling::KEEP; // Not reached.
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	if (args_info.chunk_size_arg <= 0)
	{
		std::cerr << "Chunk size must be positive." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (args_info.min_path_length_arg < 0)
	{
		std::cerr << "Minimum path length must be non-negative." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	if (args_info.reduce_samples_given && (!args_info.generated_path_count_given))
	{
		std::cerr << "'--generated-path-count' option required." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	std::size_t generated_path_count(0);
	if (args_info.generated_path_count_given)
	{
		if (args_info.generated_path_count_arg <= 0)
		{
			std::cerr << "Generated path count must be positive." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		generated_path_count = args_info.generated_path_count_arg;
	}
	
#ifndef NDEBUG
	std::cerr << "All assertions have been enabled." << std::endl;
#endif
	
	if (args_info.print_invocation_flag)
	{
		std::cerr << "Invocation:";
		for (std::size_t i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << std::endl;
	}
	
	// libdispatch on macOS does not need pthread_workqueue.
#ifdef __linux__
	pthread_workqueue_init_np();
#endif
	
	// FIXME: better alternative?
	v2m::gzip_sink_impl::init();
	
	// Window size change is going to be handled later.
	signal(SIGWINCH, SIG_IGN);

	v2m::generate_haplotypes(
		args_info.reference_arg,
		args_info.variants_arg,
		args_info.output_reference_arg,
		args_info.report_file_arg,
		args_info.null_allele_seq_arg,
		args_info.chunk_size_arg,
		args_info.min_path_length_arg,
		generated_path_count,
		sv_handling_method(args_info.structural_variants_arg),
		args_info.overwrite_flag,
		!args_info.no_check_ref_flag,
		args_info.reduce_samples_flag,
		args_info.print_subgraph_handling_flag,
		args_info.compress_output_flag
	);
		
	cmdline_parser_free(&args_info);
	
	dispatch_main();
	// Not reached b.c. pthread_exit() is eventually called in dispatch_main().
	abort();
	return EXIT_FAILURE;
}
