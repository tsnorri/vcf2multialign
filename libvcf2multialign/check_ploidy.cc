/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/check_ploidy.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void check_ploidy(lb::vcf_reader &vcf_reader, ploidy_map &out_ploidy)
	{
		out_ploidy.clear();
		vcf_reader.reset();
		vcf_reader.set_parsed_fields(lb::vcf_field::ALL);
		
		vcf_reader.fill_buffer();
		if (!vcf_reader.parse([&vcf_reader, &out_ploidy](lb::transient_variant const &var) -> bool {
			auto const *gt_field(get_variant_format(var).gt);
			for (auto const &kv : vcf_reader.sample_names())
			{
				auto const sample_no(kv.second);
				auto const &sample(var.samples()[sample_no]);
				out_ploidy[sample_no] = (*gt_field)(sample).size();
			}
		
			return false;
		}))
		{
			libbio_fail("Unable to read the first variant");
		}
	}
}
