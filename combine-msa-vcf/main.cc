/*
 * Copyright (c) 2019–2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"
#include "combine_msa.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	v2m::msa_variant_output msa_variant_output(enum enum_msa_variants const msa_variants_arg)
	{
		switch (msa_variants_arg)
		{
			case msa_variants__NULL:
			case msa_variants_arg_none:
				return v2m::msa_variant_output::NONE;
				
			case msa_variants_arg_all:
				return v2m::msa_variant_output::ALL;
				
			case msa_variants_arg_altMINUS_matchesMINUS_ref:
				return v2m::msa_variant_output::ALT_MATCHES_REF;
		}
	}
	
	
	void combine_msa(gengetopt_args_info const &args_info, std::ostream &os)
	{
		// Read the input FASTAs.
		v2m::vector_type ref_seq, alt_seq;
		v2m::read_single_fasta_seq(args_info.ref_arg, ref_seq, nullptr);
		v2m::read_single_fasta_seq(args_info.alt_arg, alt_seq, nullptr);
		
		// Combine.
		v2m::combine_msa(
			ref_seq,
			alt_seq,
			args_info.variants_arg,
			args_info.regions_arg,
			args_info.chr_arg,
			args_info.output_chr_arg,
			args_info.ploidy_arg,
			msa_variant_output(args_info.msa_variants_arg),
			os
		);
	}
}


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
	
	if (args_info.variants_arg && (!args_info.chr_arg))
	{
		std::cerr << "Chromosome identifier for input needs to be specified if variants are used for input." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	
	// Open the output file.
	if (args_info.output_given)
	{
		lb::file_ostream os;
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(false ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		lb::open_file_for_writing(args_info.output_arg, os, mode);
		combine_msa(args_info, os);
	}
	else
	{
		combine_msa(args_info, std::cout);
	}
	
	return EXIT_SUCCESS;
}
