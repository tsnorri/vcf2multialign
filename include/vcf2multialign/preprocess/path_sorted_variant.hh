/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_PATH_SORTED_VARIANT_HH
#define VCF2MULTIALIGN_PREPROCESS_PATH_SORTED_VARIANT_HH

#include <atomic>
#include <cstddef>
#include <vcf2multialign/preprocess/types.hh>


namespace vcf2multialign {

	class path_sorted_variant
	{
	protected:
		path_vector			m_paths_by_sample;
		sample_map			m_samples_by_path;
		sample_vector		m_representatives_by_path;
		
	public:
		void reserve_memory_for_paths(std::size_t path_count) { m_samples_by_path.resize(path_count); }
		void reserve_memory_for_representatives(std::size_t path_count) { m_representatives_by_path.resize(path_count); }
		
		std::size_t path_count() const { return m_samples_by_path.size(); }
		void set_paths_by_sample(path_vector const &paths) { m_paths_by_sample = paths; } // Copy.
		path_vector const &paths_by_sample() const { return m_paths_by_sample; }
		
		sample_map const &samples_by_path() const { return m_samples_by_path; }
		sample_map &samples_by_path() { return m_samples_by_path; }
		sample_vector const &representatives_by_path() const { return m_representatives_by_path; }
		sample_vector &representatives_by_path() { return m_representatives_by_path; }
		
		void invert_paths_by_sample();
		void determine_representatives_for_each_sample();
	};
}

#endif
