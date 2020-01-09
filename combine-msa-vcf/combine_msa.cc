/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include "combine_msa.hh"

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void combine_msa(
		char const *ref_path,
		char const *alt_path,
		char const *variants_path,
		char const *output_chr,
		std::uint16_t const ploidy,
		std::ostream &os
	)
	{
		// Read the input FASTAs.
		v2m::vector_type ref_seq, alt_seq;
		{
			lb::mmap_handle <char> ref_handle, alt_handle;
			ref_handle.open(ref_path);
			alt_handle.open(alt_path);
			v2m::read_single_fasta_seq(ref_handle, ref_seq, nullptr);
			v2m::read_single_fasta_seq(alt_handle, alt_seq, nullptr);
		}

		v2m::vcf_record_generator gen;
		v2m::msa_combiner combiner(output_chr, ploidy, os);

		if (variants_path)
		{
			gen.open_variants_file(variants_path);
			gen.prepare();
			gen.vcf_reader().set_parsed_fields(lb::vcf_field::ALL);
		}

		combiner.handle_msa(ref_seq, alt_seq, gen);
	}
}
