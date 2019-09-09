/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_SAMPLE_INDEXER_HH
#define VCF2MULTIALIGN_PREPROCESS_SAMPLE_INDEXER_HH

#include <cstddef>
#include <cstdint>
#include <tuple>


namespace vcf2multialign {

	// Map sample numbers and chromosome indices to consecutive integers.
	// Assume same ploidy for each donor.
	class sample_indexer
	{
	public:
		typedef std::size_t								sample_type;
		typedef	std::size_t								donor_type;
		typedef std::uint8_t							chr_type;
		typedef std::tuple <donor_type, chr_type>		donor_and_chr_pair;
		
	protected:
		donor_type	m_donor_count{};
		chr_type	m_chr_count{};

	public:
		sample_indexer() = default;
		
		sample_indexer(donor_type const donor_count, chr_type const chr_count):
			m_donor_count(donor_count),
			m_chr_count(chr_count)
		{
		}
		
		inline void set_counts(donor_type const donor_count, chr_type const chr_count) { m_donor_count = donor_count; m_chr_count = chr_count; }
		inline donor_type donor_count() const { return m_donor_count; }
		inline chr_type chr_count() const { return m_chr_count; }
		inline sample_type total_samples() const { return m_donor_count * m_chr_count; }
		
		inline sample_type sample_idx(donor_type const donor_idx, chr_type const chr_idx) const
		{
			return donor_idx * m_chr_count + chr_idx;
		}
		
		inline donor_and_chr_pair donor_and_chr_idx(std::size_t sample_idx) const
		{
			return donor_and_chr_pair(sample_idx / m_chr_count, sample_idx % m_chr_count);
		}
	};
}

#endif
