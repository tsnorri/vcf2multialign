/*
 * Copyright (c) 2019–2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"
#include "combine_msa.hh"


namespace v2m	= vcf2multialign;


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	if (args_info.show_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}
	
	if (args_info.ploidy_arg <= 0)
	{
		std::cerr << "Ploidy needs to be positive." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	
	// Read the input FASTAs.
	v2m::vector_type ref_seq, alt_seq;
	v2m::read_single_fasta_seq(args_info.ref_arg, ref_seq, nullptr);
	v2m::read_single_fasta_seq(args_info.alt_arg, alt_seq, nullptr);
	
	// Combine.
	v2m::combine_msa(
		ref_seq,
		alt_seq,
		args_info.variants_arg,
		args_info.output_chr_arg,
		args_info.ploidy_arg,
		std::cout
	);
	
	return EXIT_SUCCESS;
}
