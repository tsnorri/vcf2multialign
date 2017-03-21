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
	
	
	void variant_base::set_gt(std::size_t const alt, std::size_t const sample_no, std::size_t const idx, bool const is_phased)
	{
		if (! (sample_no < m_samples.size()))
			m_samples.resize(sample_no + 1);
		
		auto &sample(m_samples[sample_no]);
		
		if (! (idx < sample.genotype.size()))
			sample.genotype.resize(idx + 1);
		
		auto &gt(sample.genotype[idx]);
		
		gt.alt = alt;
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
