/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <range/v3/all.hpp>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>


namespace vcf2multialign {
	
	void path_sorted_variant::invert_paths_by_sample()
	{
		for (auto const &tup : ranges::view::enumerate(m_paths_by_sample))
		{
			auto const [sample_idx, path_idx] = tup;
			m_samples_by_path[path_idx].push_back(sample_idx);
		}
	}
}
