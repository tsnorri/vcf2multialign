/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_SAMPLE_SORTER_HH
#define VCF2MULTIALIGN_PREPROCESS_SAMPLE_SORTER_HH

#include <libbio/packed_vector.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/vcf_subfield.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vector>


namespace vcf2multialign {
	
	class sample_sorter
	{
	protected:
		typedef std::uint16_t								path_number_type;
		typedef libbio::packed_vector <2>					branch_vector;

	protected:
		sample_indexer const		*m_sample_indexer{};
		libbio::vcf_info_field_end	*m_end_field{};
		path_vector					m_src_paths;
		path_vector					m_dst_paths;
		path_vector					m_branching_paths;			// Path numbers for branching paths.
		branch_vector				m_branches_by_path_index;
		std::vector <std::size_t>	m_end_positions_by_sample;
		path_number_type			m_path_counter{};

	public:
		sample_sorter() = default;
	
		sample_sorter(libbio::vcf_reader &reader, sample_indexer const &indexer):
			m_sample_indexer(&indexer),
			m_src_paths(indexer.total_samples(), 0),
			m_dst_paths(indexer.total_samples(), 0),
			m_branching_paths(indexer.total_samples(), 0),
			m_branches_by_path_index(indexer.total_samples()),
			m_end_positions_by_sample(indexer.total_samples())
		{
			reader.get_info_field_ptr("END", m_end_field);
		}
	
		void prepare_for_next_subgraph();
		void sort_by_variant_and_alt(libbio::variant const &var, std::uint8_t const expected_alt_idx);
		path_vector const &paths_by_sample() const { return m_src_paths; }
		path_number_type path_count() const { return m_path_counter; }
		
	protected:
		inline std::size_t variant_end_pos(libbio::variant const &var) const { return libbio::variant_end_pos(var, *m_end_field); }
		
		bool check_variant_for_sample_and_update_state(
			std::size_t const sample_idx,
			libbio::variant const &var,
			std::uint8_t const alt_idx
		);
	};
}

#endif
