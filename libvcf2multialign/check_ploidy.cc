/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/check_ploidy.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;
namespace vcf	= libbio::vcf;


namespace vcf2multialign {
	
	void check_ploidy(vcf::reader &vcf_reader, ploidy_map &out_ploidy)
	{
		out_ploidy.clear();
		vcf_reader.set_parsed_fields(vcf::field::ALL);
		
		bool handled_variant(false);
		vcf_reader.parse([&vcf_reader, &out_ploidy, &handled_variant](vcf::transient_variant const &var) -> bool {
			handled_variant = true;
			auto const *gt_field(get_variant_format(var).gt);
			for (auto const &kv : vcf_reader.sample_indices_by_name()) // FIXME: use sample_names_by_index instead.
			{
				auto const sample_no(kv.second);
				libbio_assert_lt(0, sample_no);
				libbio_assert_lte(sample_no, var.samples().size());
				auto const &sample(var.samples()[sample_no - 1]);
				out_ploidy[sample_no] = (*gt_field)(sample).size();
			}
		
			return false;
		});
		
		if (!handled_variant)
		{
			libbio_fail("Unable to read the first variant");
		}
	}
}
