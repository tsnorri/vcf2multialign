/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/util.hh>
#include <vcf2multialign/variant.hh>


namespace vcf2multialign {

	size_t variant_base::zero_based_pos() const
	{
		always_assert(0 != m_pos, "Unexpected position");
		return m_pos - 1;
	}
	
	
	void variant_base::set_gt(std::size_t const alt, std::size_t const sample, std::size_t const idx, bool const is_phased)
	{
		if (! (sample < m_gt.size()))
			m_gt.resize(sample + 1);
		
		auto gt(m_gt[sample]);
		if (! (idx < gt.alts.size()))
			gt.alts.resize(idx + 1);
		
		gt.alts[idx] = alt;
		gt.is_phased = is_phased;
	}
	
	
	void transient_variant::reset()
	{
		std::string_view empty(nullptr, 0);
		m_chrom_id = empty;
		m_ref = empty;
	}
	
	
	void variant::reset()
	{
		m_chrom_id.clear();
		m_ref.clear();
	}
}
