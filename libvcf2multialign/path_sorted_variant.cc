/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <range/v3/all.hpp>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>


namespace vcf2multialign {
	
	void path_sorted_variant::invert_paths_by_sample()
	{
		for (auto const &[sample_idx, path_idx] : ranges::view::enumerate(m_paths_by_sample))
			m_samples_by_path[path_idx].push_back(sample_idx);
	}
	
	void path_sorted_variant::determine_representatives_for_each_sample()
	{
		std::fill(m_representatives_by_path.begin(), m_representatives_by_path.end(), SAMPLE_NUMBER_MAX);
		for (auto const &[sample_idx, path_idx] : ranges::view::enumerate(m_paths_by_sample))
		{
			if (SAMPLE_NUMBER_MAX == m_representatives_by_path[path_idx])
				m_representatives_by_path[path_idx] = sample_idx;
		}
	}
}
