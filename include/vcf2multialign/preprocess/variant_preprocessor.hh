/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_VARIANT_PREPROCESSOR_HH
#define VCF2MULTIALIGN_PREPROCESS_VARIANT_PREPROCESSOR_HH

#include <libbio/matrix.hh>
#include <libbio/packed_matrix.hh>
#include <libbio/vcf/variant.hh>
#include <ostream>
#include <string>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/sample_sorter.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vcf2multialign/types.hh>
#include <vector>


namespace vcf2multialign {
	
	class variant_preprocessor
	{
	protected:
		typedef std::vector <libbio::variant>	variant_stack;
	
	protected:
		libbio::vcf_reader						*m_reader{};
		std::ostream							*m_ostream{};
		vector_type const						*m_reference{};
		libbio::vcf_info_field_end				*m_end_field{};
		std::string								m_chromosome_name;
		variant_stack							m_overlapping_variants;
		sample_indexer							m_sample_indexer;
		sample_sorter							m_sample_sorter;
		libbio::matrix <std::uint8_t>			m_output_gt;
		std::vector <std::string>				m_sample_names;
		std::size_t								m_overlap_start{};
		std::size_t								m_overlap_end{};
		std::size_t								m_minimum_subgraph_distance{};
		std::size_t								m_output_lineno{};
		
	public:
		variant_preprocessor(
			libbio::vcf_reader &reader,
			std::ostream &ostream,
			vector_type const &reference,
			std::string const chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count
		):
			m_reader(&reader),
			m_ostream(&ostream),
			m_reference(&reference),
			m_chromosome_name(chr_name),
			m_sample_indexer(donor_count, chr_count),
			m_sample_sorter(reader, m_sample_indexer),
			m_output_gt(chr_count, donor_count)
		{
			reader.get_info_field_ptr("END", m_end_field);
		}
		
		void process_and_output(std::vector <std::string> const &field_names_for_filter_by_assigned);
		
	protected:
		void update_sample_names();
		void output_headers();
		
		void process_overlap_stack();
		void output_sorted(path_sorted_variant const &psv, libbio::packed_matrix <1> const &alt_sequences_by_path);
		
		void output_single_variant(libbio::variant const &var);
		void output_combined_variants(
			path_sorted_variant const &psv,
			std::vector <std::string> const &alt_strings_by_path
		);
		void output_sorted_handle_paths(
			std::size_t const row_idx,
			libbio::variant const &var,
			std::size_t const next_var_pos,
			std::string_view const &reference_sv,
			libbio::packed_matrix <1> const &alt_sequences_by_path,
			sample_map const &samples_by_path,
			std::vector <std::string> &alt_strings_by_path,
			std::vector <std::size_t> &alt_string_output_positions
		) const;
	};
}

#endif
