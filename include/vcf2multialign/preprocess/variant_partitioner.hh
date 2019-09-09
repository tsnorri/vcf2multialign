/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_VARIANT_PARTITIONER_HH
#define VCF2MULTIALIGN_PREPROCESS_VARIANT_PARTITIONER_HH

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
#include <vcf2multialign/variant_processor_delegate.hh>
#include <vector>


namespace vcf2multialign {
	class variant_partitioner; // Fwd.
}


namespace vcf2multialign {
	
	struct cut_position_list
	{
		std::vector <std::size_t>	handled_line_numbers;	// In addition to the cut positions, store the line numbers that were handled.
		position_vector				positions;
		path_number_type			max_segment_size{};
		std::size_t					donor_count{};
		std::uint8_t				chr_count{};
		
		// Ignore the version for now.
		template <typename t_archive>
		void serialize(t_archive &archive, std::uint32_t const version)
		{
			archive(handled_line_numbers, positions, max_segment_size, donor_count, chr_count);
		}
	};
	
	
	struct variant_partitioner_delegate : public variant_processor_delegate, public sample_sorter_delegate {};
	
	
	class variant_partitioner
	{
	protected:
		struct dp_ctx; // Fwd
		
	protected:
		variant_partitioner_delegate					*m_delegate{};
		libbio::vcf_reader								*m_reader{};
		libbio::vcf_info_field_end						*m_end_field{};
		vector_type const								*m_reference{};
		sample_indexer									m_sample_indexer;
		std::string										m_chromosome_name;
		std::size_t										m_minimum_subgraph_distance{};
		std::atomic_size_t								m_processed_count{};
		
	public:
		variant_partitioner(
			variant_partitioner_delegate &delegate,
			libbio::vcf_reader &reader,
			vector_type const &reference,
			std::string const chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count,
			std::size_t const minimum_subgraph_distance
		):
			m_delegate(&delegate),
			m_reader(&reader),
			m_end_field(reader.get_end_field_ptr()),
			m_reference(&reference),
			m_sample_indexer(donor_count, chr_count),
			m_chromosome_name(chr_name),
			m_minimum_subgraph_distance(minimum_subgraph_distance)
		{
		}
		
		bool partition(
			std::vector <std::string> const &field_names_for_filter_by_assigned,
			cut_position_list &out_cut_positions
		);
		
		std::size_t processed_count() const { return m_processed_count.load(std::memory_order_relaxed); }
		
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
			std::size_t		previous_idx{SIZE_MAX};	// in cut position tree.
			
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
			
			dp_ctx(sample_sorter_delegate &delegate, libbio::vcf_reader &reader, sample_indexer const &indexer):
				sorter(delegate, reader, indexer)
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
