/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <numeric>
#include <range/v3/all.hpp>
#include "variant_writer.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void variant_writer::output_vcf_header() const
	{
		auto &os(*m_os);
		os << "##fileformat=VCFv4.3\n";
		os << "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n";
		os << "##FILTER=<ID=ALT_EQ_TO_REF,Description=\"Variant called by the VC is equivalent to the reference\">\n";
		os << "##FILTER=<ID=GT_NOT_SET,Description=\"All GT values are equal to zero\">\n";
		os << "##INFO=<ID=OC,Number=1,Type=Integer,Description=\"Number of overlapping VC variants\">\n";
		os << "##INFO=<ID=USRA,Number=0,Type=Flag,Description=\"Uses source reference for ALT (ALT matched REF after projection and was substituted with a substring of the original reference)\">\n";
		
		if (m_should_output_msa_deduced_variants)
			os << "##INFO=<ID=VS,Number=1,Type=String,Description=\"Variant source (MSA for multiple sequence alignment, VC for variant caller)\">\n";
		
		os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
	}
	
	
	void variant_writer::handle_variant_description(variant_description &&desc)
	{
		if (!m_should_output_msa_deduced_variants && variant_origin::MSA == desc.origin)
			return;
		
		auto &os(*m_os);
		
		// CHROM
		os << m_output_chr_id << '\t';
			
		// POS, ID, REF
		libbio_assert_lt(0, desc.ref.size());
		os << (1 + desc.position) << "\t.\t" << desc.ref << '\t';
			
		// ALT
		if (desc.alt.empty())
			os << "<DEL>";
		else
			os << desc.alt;
			
		// QUAL
		os << "\t.\t";
			
		// FILTER
		if (desc.filters.empty())
			os << "PASS";
		else
			ranges::copy(desc.filters, ranges::make_ostream_joiner(os, ";"));
			
		// INFO
		os << "\tOC=" << desc.overlap_count;
		if (m_should_output_msa_deduced_variants)
		{
			os << ";VS=";
			switch (desc.origin)
			{
				case variant_origin::MSA:
				{
					os << "MSA";
					break;
				}
					
				case variant_origin::VC:
				{
					os << "VC";
					break;
				}
			}
		}
		
		if (desc.had_alt_eq_to_ref)
			os << ";USRA";
		
		// FORMAT
		os << "\tGT\t";
			
		// Sample.
		ranges::copy(desc.genotype | rsv::transform([](auto const gt) -> std::size_t { return gt; }), ranges::make_ostream_joiner(os, "/"));
			
		os << '\n';
	}
}
