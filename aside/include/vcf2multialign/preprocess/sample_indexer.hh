/*
 * Copyright (c) 2019-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_SAMPLE_INDEXER_HH
#define VCF2MULTIALIGN_PREPROCESS_SAMPLE_INDEXER_HH

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vcf2multialign/variant_format.hh>


namespace vcf2multialign {

	// Map given sample numbers and chromosome indices to consecutive integers by using one variant as input.
	class sample_indexer
	{
	public:
		typedef std::size_t								sample_type;
		typedef	std::size_t								donor_type;
		typedef std::uint8_t							chr_type;
		typedef std::tuple <donor_type, chr_type>		donor_and_chr_pair;
		
	protected:
		std::vector <chr_type>	m_cumulative_sum{0, 1};

	public:
		template <typename t_variant>
		void prepare(t_variant const &variant)
		{
			auto const &samples(variant.samples());
			m_cumulative_sum.clear();
			m_cumulative_sum.resize(1 + samples.size(), 0);
			
			auto const *gt_field(get_variant_format(variant).gt);
			libbio_assert(gt_field);
			for (auto const &[idx, sample] : ranges::views::enumerate(samples))
			{
				std::vector <libbio::vcf::sample_genotype> const &gt((*gt_field)(sample));
				m_cumulative_sum[1 + idx] = m_cumulative_sum[idx] + gt.size();
			}
		}
		
		sample_indexer() = default;
		
		template <typename t_variant>
		sample_indexer(t_variant const &variant)
		{
			// Try to determine that the variant is actually valid.
			libbio_assert_lt_msg(0, variant.lineno(), "Expected the passed variant to be valid");
			prepare(variant);
		}
		
		inline sample_type total_samples() const { return m_cumulative_sum.back(); }
		
		inline sample_type sample_idx(donor_type const donor_idx, chr_type const chr_idx) const
		{
			libbio_assert_lte(0, donor_idx);
			libbio_assert_lt(donor_idx + 1, m_cumulative_sum.size());
			auto const retval(m_cumulative_sum[donor_idx] + chr_idx);
			libbio_assert_lt(retval, m_cumulative_sum[1 + donor_idx]);
			return retval;
		}
		
		inline donor_and_chr_pair donor_and_chr_idx(std::size_t sample_idx) const
		{
			auto const begin(m_cumulative_sum.begin());
			auto const end(m_cumulative_sum.end());
			auto const it(std::upper_bound(begin, end, sample_idx));
			libbio_assert_neq(it, end); // Since we have the cumulative sum of the counts, the index should always be less than the maximum value.
			libbio_assert_neq(it, begin);
			return {std::distance(begin, it) - 1, sample_idx - *(it - 1)};
		}
	};
}

#endif
