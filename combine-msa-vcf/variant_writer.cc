/*
 * Copyright (c) 2020 Tuukka Norri
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
	
	void variant_writer::merge_output_variants(std::size_t const partition_point)
	{
		// Merge the partitions of sorted variants.
		std::inplace_merge(
			m_output_variants.begin(),
			m_output_variants.begin() + partition_point,
			m_output_variants.end(),
			[](auto const &lhs, auto const &rhs){
				return lhs.position < rhs.position;
			}
		);
		libbio_assert(
			std::is_sorted(
				m_output_variants.begin(),
				m_output_variants.end(),
				[](auto const &lhs, auto const &rhs){ return lhs.position < rhs.position; }
			)
		);
	}
	
	
	void variant_writer::output_vcf_header() const
	{
		auto &os(*m_os);
		os << "##fileformat=VCFv4.3\n";
		os << "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n";
		os << "##FILTER=<ID=ALT_EQ_TO_REF,Description=\"Variant called by the VC is equivalent to the reference\">\n";
		os << "##FILTER=<ID=GT_NOT_SET,Description=\"All GT values are equal to zero.\">\n";
		os << "##INFO=<ID=OC,Number=1,Type=Integer,Description=\"Number of overlapping VC variants\">\n";
		os << "##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Variant source (MSA for multiple sequence alignment, VC for variant caller)\">\n";
		os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
	}
	
	
	void variant_writer::filter_processed_variants_and_output(std::size_t const min_unhandled_ref_pos)
	{
		// Omit filtering for now, except for checking REF against ALT.
		libbio_assert(
			std::is_sorted(m_output_variants.begin(), m_output_variants.end(), [](auto const &lhs, auto const &rhs){
				return lhs.position < rhs.position;
			})
		);
		auto &os(*m_os);
		std::vector <std::string> filters;
		auto var_it(m_output_variants.begin());
		auto const var_end(m_output_variants.end());
		while (var_it != var_end)
		{
			auto const &desc(*var_it);
			if (min_unhandled_ref_pos <= desc.position)
				break;
			
			filters.clear();
			if (desc.ref == desc.alt)
				filters.emplace_back("ALT_EQ_TO_REF");
			
			if (0 == std::accumulate(desc.genotype.begin(), desc.genotype.end(), std::uint16_t(0)))
				filters.emplace_back("GT_NOT_SET");
			
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
			if (filters.empty())
				os << "PASS";
			else
				ranges::copy(filters, ranges::make_ostream_joiner(os, ";"));
			
			// INFO
			os << "\tOC=" << desc.overlap_count;
			os << ";ORIGIN=";
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
			
			// FORMAT
			os << "\tGT\t";
			
			// Sample.
			ranges::copy(desc.genotype | rsv::transform([](auto const gt) -> std::size_t { return gt; }), ranges::make_ostream_joiner(os, "/"));
			
			os << '\n';
			++var_it;
		}
		
		m_output_variants.erase(m_output_variants.begin(), var_it);
	}
}
