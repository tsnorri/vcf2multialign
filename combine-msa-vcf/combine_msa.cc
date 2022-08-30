/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include "combine_msa.hh"

namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void combine_msa(
		vector_type const &ref_seq,
		vector_type const &alt_seq,
		char const *variants_path,
		char const *regions_bed_path,
		char const *input_chr,
		char const *output_chr,
		std::uint16_t const ploidy,
		std::ostream &os,
		bool const log_status
	)
	{
		v2m::vcf_record_generator gen(input_chr);
		v2m::msa_combiner combiner(output_chr, ploidy, os, log_status);

		if (variants_path)
		{
			gen.open_variants_file(variants_path, regions_bed_path);
			gen.prepare();
			gen.vcf_reader().set_parsed_fields(vcf::field::ALL);
		}
		else
		{
			gen.setup_empty_input();
		}

		combiner.process_msa(ref_seq, alt_seq, gen);
	}
}
