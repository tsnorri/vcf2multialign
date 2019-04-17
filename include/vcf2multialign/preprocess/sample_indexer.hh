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
	protected:
		std::size_t		m_donor_count{};
		std::uint8_t	m_chr_count{};
		
	public:
		sample_indexer() = default;
		
		sample_indexer(std::size_t donor_count, std::uint8_t chr_count):
			m_donor_count(donor_count),
			m_chr_count(chr_count)
		{
		}
		
		inline void set_counts(std::size_t donor_count, std::uint8_t chr_count) { m_donor_count = donor_count; m_chr_count = chr_count; }
		inline std::size_t donor_count() const { return m_donor_count; }
		inline std::uint8_t chr_count() const { return m_chr_count; }
		inline std::size_t total_samples() const { return m_donor_count * m_chr_count; }
		
		inline std::size_t sample_idx(std::size_t donor_idx, std::uint8_t chr_idx) const
		{
			return donor_idx * m_chr_count + chr_idx;
		}
		
		inline std::tuple <std::size_t, std::uint8_t> donor_and_chr_idx(std::size_t sample_idx) const
		{
			return std::make_tuple(sample_idx / m_chr_count, sample_idx % m_chr_count);
		}
	};
}

#endif
