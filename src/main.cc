/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <vcf2multialign/generate_haplotypes.hh>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif

#include "cmdline.h"


namespace v2m	= vcf2multialign;


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	if (args_info.chunk_size_arg <= 0)
	{
		std::cerr << "Chunk size must be positive." << std::endl;
		exit(EXIT_FAILURE);
	}

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
		args_info.overwrite_flag,
		!args_info.no_check_ref_flag
	);
	
	// Not reached b.c. pthread_exit() is eventually called in generate_haplotypes().
	return EXIT_SUCCESS;
}
