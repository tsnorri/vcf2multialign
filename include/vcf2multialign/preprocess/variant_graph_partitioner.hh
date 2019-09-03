/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_VARIANT_GRAPH_PARTITIONER_HH
#define VCF2MULTIALIGN_PREPROCESS_VARIANT_GRAPH_PARTITIONER_HH

#include <libbio/matrix.hh>
#include <libbio/int_matrix.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/vcf_subfield_def.hh>
#include <list>
#include <ostream>
#include <string>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/sample_sorter.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vcf2multialign/types.hh>
#include <vector>


namespace vcf2multialign {
	class variant_graph_partitioner; // Fwd.
}


namespace vcf2multialign { namespace detail {
}}


namespace vcf2multialign {

	class variant_graph_partitioner
	{
	protected:
		struct dp_ctx; // Fwd
		
	public:
		struct cut_position_list; // Fwd
		typedef std::size_t								position_type;
		typedef std::vector <position_type>				position_vector;
		
	protected:
		libbio::vcf_reader								*m_reader{};
		libbio::vcf_info_field_end						*m_end_field{};
		vector_type const								*m_reference{};
		sample_indexer									m_sample_indexer;
		std::string										m_chromosome_name;
		std::size_t										m_minimum_subgraph_distance{};
		std::size_t										m_processed_count{};
		
	protected:
		static libbio::vcf_info_field_end *get_end_field_ptr(libbio::vcf_reader &reader)
		{
			libbio::vcf_info_field_end *retval{};
			reader.get_info_field_ptr("END", retval);
			return retval;
		}

		variant_graph_partitioner(
			libbio::vcf_reader &reader,
			libbio::vcf_info_field_end &end_field,
			vector_type const &reference,
			std::string const chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count,
			std::size_t const minimum_subgraph_distance
		):
			m_reader(&reader),
			m_end_field(&end_field),
			m_reference(&reference),
			m_sample_indexer(donor_count, chr_count),
			m_chromosome_name(chr_name),
			m_minimum_subgraph_distance(minimum_subgraph_distance)
		{
		}

	public:
		variant_graph_partitioner(
			libbio::vcf_reader &reader,
			vector_type const &reference,
			std::string const chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count,
			std::size_t const minimum_subgraph_distance
		):
			variant_graph_partitioner(
				reader,
				*variant_graph_partitioner::get_end_field_ptr(reader),
				reference,
				chr_name,
				donor_count,
				chr_count,
				minimum_subgraph_distance
			)
		{
		}
		
		bool partition(
			std::vector <std::string> const &field_names_for_filter_by_assigned,
			cut_position_list &out_cut_positions
		);
		
	public:
		struct cut_position_list
		{
			position_vector		positions;
			path_number_type	max_segment_size{};
			
			// Ignore the version for now.
			template <typename t_archive>
			void serialize(t_archive &archive, std::uint32_t const version) { archive(positions, max_segment_size); }
		};
		
	protected:
		struct cut_position; // Fwd
		
		void check_closable(
			std::size_t const var_pos,
			std::vector <cut_position> const &cut_position_tree,
			std::list <dp_ctx> &unclosable_partitions,
			std::list <dp_ctx> &closable_partitions
		);
		
		struct cut_position
		{
			position_type	value{};
			std::size_t		previous_idx{SIZE_MAX};
			
			cut_position() = default;
			
			cut_position(position_type const value_, std::size_t const previous_idx_):
				value(value_),
				previous_idx(previous_idx_)
			{
			}
		};
		
		struct dp_ctx
		{
			sample_sorter		sorter;
			std::size_t			start_position_idx{};
			path_number_type	max_size{};
			
			dp_ctx() = default;
			
			dp_ctx(libbio::vcf_reader &reader, sample_indexer const &indexer):
				sorter(reader, indexer)
			{
				sorter.prepare_for_next_subgraph();
			}
		
			// Called once per variant only.
			void chain_previous(dp_ctx const &other, position_type const pos, std::vector <cut_position> &cut_position_tree)
			{
				cut_position_tree.emplace_back(pos, other.start_position_idx);
				start_position_idx = cut_position_tree.size() - 1;
				max_size = other.max_size;
			}
			
			void count_paths(libbio::transient_variant const &var, std::size_t const alt_count)
			{
				for (std::size_t i(0); i < alt_count; ++i)
				{
					sorter.sort_by_variant_and_alt(var, 1 + i);
					max_size = std::max(max_size, sorter.path_count());
				}
			}

			void output_reversed_path(std::ostream &, std::vector <cut_position> const &) const;
		};
		
		friend std::ostream &operator<<(std::ostream &, dp_ctx const &);
	};
}

#endif
