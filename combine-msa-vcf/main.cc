/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "cmdline.h"
#include "msa_combiner.hh"
#include "vcf_record_generator.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
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
	{
		lb::mmap_handle <char> ref_handle, alt_handle;
		ref_handle.open(args_info.ref_arg);
		alt_handle.open(args_info.alt_arg);
		v2m::read_single_fasta_seq(ref_handle, ref_seq, nullptr);
		v2m::read_single_fasta_seq(alt_handle, alt_seq, nullptr);
	}
	
	v2m::vcf_record_generator gen;
	v2m::msa_combiner combiner(args_info.output_chr_arg, args_info.ploidy_arg);
	
	gen.open_variants_file(args_info.variants_arg);
	gen.prepare();
	gen.vcf_reader().set_parsed_fields(lb::vcf_field::ALL);
	
	combiner.handle_msa(ref_seq, alt_seq, gen);
	
	return EXIT_SUCCESS;
}
