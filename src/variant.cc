/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/util.hh>
#include <vcf2multialign/variant.hh>


namespace vcf2multialign {

	std::size_t variant_base::zero_based_pos() const
	{
		always_assert(0 != m_pos, "Unexpected position");
		return m_pos - 1;
	}
	
	
	void variant_base::set_gt(std::size_t const alt, std::size_t const sample_no, std::size_t const idx, bool const is_phased)
	{
		// Check that the samples are given in consecutive order.
		always_assert(0 != sample_no);
		always_assert(sample_no == 1 || sample_no == m_sample_count - 1 || sample_no == m_sample_count);
		m_sample_count = 1 + sample_no;
		
		if (! (m_sample_count <= m_samples.size()))
			m_samples.resize(m_sample_count);
		
		auto &sample(m_samples[sample_no]);
		
		// Again check the order.
		always_assert(0 == idx || idx == sample.m_gt_count);
		sample.m_gt_count = 1 + idx;
		
		if (! (sample.m_gt_count <= sample.m_genotype.size()))
			sample.m_genotype.resize(sample.m_gt_count);
		
		auto &gt(sample.m_genotype[idx]);
		
		gt.alt = alt;
		gt.is_phased = is_phased;
	}
	
	
	void variant_base::set_alt_sv_type(sv_type const svt, std::size_t const pos)
	{
		if (! (pos < m_alt_sv_types.size()))
			m_alt_sv_types.resize(pos + 1);
		
		m_alt_sv_types[pos] = svt;
	}
	
	
	void transient_variant::reset()
	{
		superclass::reset();
		std::string_view empty(nullptr, 0);
		m_chrom_id = empty;
		m_ref = empty;
	}
	
	
	void variant::reset()
	{
		superclass::reset();
		m_chrom_id.clear();
		m_ref.clear();
	}
	
	
	std::ostream &operator<<(std::ostream &os, sample_field const &sample_field)
	{
		bool is_first(true);
		for (auto const &gt_field : sample_field.get_genotype())
		{
			if (!is_first)
				os << (gt_field.is_phased ? '|' : '/');
			
			is_first = false;
			os << gt_field.alt;
		}
		
		return os;
	}
}
