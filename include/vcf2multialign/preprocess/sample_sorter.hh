/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_SAMPLE_SORTER_HH
#define VCF2MULTIALIGN_PREPROCESS_SAMPLE_SORTER_HH

#include <libbio/int_vector.hh>
#include <libbio/vcf/subfield.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/variant_end_pos.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vector>


namespace vcf2multialign {
	
	struct sample_sorter_delegate
	{
		virtual ~sample_sorter_delegate() {}
		virtual void sample_sorter_found_overlapping_variant(libbio::vcf::variant const &var, std::size_t const sample_idx, std::size_t const prev_end_pos) = 0;
	};
	
	
	class sample_sorter
	{
	protected:
		typedef std::uint16_t								path_number_type;
		typedef libbio::atomic_int_vector <2>				branch_vector;

	protected:
		sample_sorter_delegate		*m_delegate{};
		sample_indexer const		*m_sample_indexer{};
		libbio::vcf::info_field_end	*m_end_field{};
		path_vector					m_src_paths;
		path_vector					m_dst_paths;
		path_vector					m_branching_paths;			// Path numbers for branching paths.
		branch_vector				m_branches_by_path_index;
		std::vector <std::size_t>	m_end_positions_by_sample;
		path_number_type			m_path_counter{};

	public:
		sample_sorter() = default;
	
		sample_sorter(sample_sorter_delegate &delegate, libbio::vcf::reader &reader, sample_indexer const &indexer):
			m_delegate(&delegate),
			m_sample_indexer(&indexer),
			m_src_paths(indexer.total_samples(), 0),
			m_dst_paths(indexer.total_samples(), 0),
			m_branching_paths(indexer.total_samples(), 0),
			m_branches_by_path_index(indexer.total_samples()),
			m_end_positions_by_sample(indexer.total_samples())
		{
			reader.get_info_field_ptr("END", m_end_field);
		}
	
		void set_delegate(sample_sorter_delegate &delegate) { m_delegate = &delegate; }
		void set_sample_indexer(sample_indexer &indexer) { m_sample_indexer = &indexer; }
		void prepare_for_next_subgraph();
		void sort_by_variant_and_alt(libbio::vcf::variant const &var, std::uint8_t const expected_alt_idx);
		path_vector const &paths_by_sample() const { return m_src_paths; }
		path_number_type path_count() const { return m_path_counter; }
		
	protected:
		inline std::size_t variant_end_pos(libbio::vcf::variant const &var) const { return libbio::vcf::variant_end_pos(var, *m_end_field); }
		
		bool check_variant_for_sample_and_update_state(
			std::size_t const sample_idx,
			libbio::vcf::variant const &var,
			std::uint8_t const alt_idx
		);
	};
}

#endif
