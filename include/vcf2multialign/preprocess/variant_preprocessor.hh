/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_VARIANT_PREPROCESSOR_HH
#define VCF2MULTIALIGN_PREPROCESS_VARIANT_PREPROCESSOR_HH

#include <libbio/matrix.hh>
#include <libbio/int_matrix.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/vcf_subfield_def.hh>
#include <ostream>
#include <queue>
#include <string>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/sample_sorter.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vcf2multialign/preprocess/variant_graph.hh>
#include <vcf2multialign/types.hh>
#include <vector>


namespace vcf2multialign {
	class variant_preprocessor; // Fwd.
}


namespace vcf2multialign { namespace detail {
	
	struct overlap_stack_entry
	{
		libbio::variant const	*variant{};
		std::size_t				node_number{};
		std::size_t				max_alt_edge_aligned_dst_pos{};
		
		overlap_stack_entry(libbio::variant const &variant_, std::size_t node_number_, std::size_t max_alt_edge_aligned_dst_pos_):
			variant(&variant_),
			node_number(node_number_),
			max_alt_edge_aligned_dst_pos(max_alt_edge_aligned_dst_pos_)
		{
		}
	};
	
	// Use with a priority queue to sort the entries by variant end position in ascending order.
	struct overlap_stack_compare
	{
		libbio::vcf_info_field_end const *end_field{};
		
		overlap_stack_compare(libbio::vcf_info_field_end const &end_field_):
			end_field(&end_field_)
		{
		}
		
		bool operator()(overlap_stack_entry const &lhs, overlap_stack_entry const &rhs)
		{
			// priority_queue sorts in descending order by default, hence the use of operator >.
			return libbio::variant_end_pos(*lhs.variant, *end_field) > libbio::variant_end_pos(*rhs.variant, *end_field);
		}
	};
}}


namespace vcf2multialign {
	
	struct variant_preprocessor_delegate
	{
		virtual ~variant_preprocessor_delegate() {}
		virtual void variant_preprocessor_no_field_for_identifier(std::string const &identifier) = 0;
		virtual void variant_preprocessor_found_variant_with_position_greater_than_reference_length(libbio::transient_variant const &var) = 0;
		virtual void variant_preprocessor_found_variant_with_no_suitable_alts(libbio::transient_variant const &var) = 0;
		virtual void variant_preprocessor_found_filtered_variant(libbio::transient_variant const &var, libbio::vcf_info_field_base const &field) = 0;
		virtual void variant_preprocessor_found_variant_with_ref_mismatch(libbio::transient_variant const &var, std::string_view const &ref_sub) = 0;
		virtual void variant_preprocessor_will_handle_subgraph(std::size_t const variant_count, std::size_t const path_count) = 0;
	};
	
	
	class variant_preprocessor
	{
	protected:
		typedef std::vector <libbio::variant>			variant_stack;
		typedef std::priority_queue <
			detail::overlap_stack_entry,
			std::vector <detail::overlap_stack_entry>,
			detail::overlap_stack_compare
		>												overlap_stack;
	
	protected:
		libbio::vcf_reader								*m_reader{};
		vector_type const								*m_reference{};
		libbio::vcf_info_field_end						*m_end_field{};
		variant_preprocessor_delegate					*m_delegate{};
		variant_graph									m_graph;
		std::string										m_chromosome_name;
		variant_stack									m_subgraph_variants;
		sample_indexer									m_sample_indexer;
		sample_sorter									m_sample_sorter;
		std::vector <std::string>						m_sample_names;
		overlap_stack									m_overlap_stack;
		std::vector <std::size_t>						m_unhandled_alt_csum;
		std::size_t										m_minimum_subgraph_distance{};
		std::size_t										m_output_lineno{};
		std::size_t										m_processed_count{};
	
	protected:
		static libbio::vcf_info_field_end *get_end_field_ptr(libbio::vcf_reader &reader)
		{
			libbio::vcf_info_field_end *retval{};
			reader.get_info_field_ptr("END", retval);
			return retval;
		}
		
		variant_preprocessor(
			variant_preprocessor_delegate &delegate,
			libbio::vcf_reader &reader,
			libbio::vcf_info_field_end &end_field,
			vector_type const &reference,
			std::string const chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count
		):
			m_reader(&reader),
			m_reference(&reference),
			m_end_field(&end_field),
			m_delegate(&delegate),
			m_chromosome_name(chr_name),
			m_sample_indexer(donor_count, chr_count),
			m_sample_sorter(reader, m_sample_indexer),
			m_overlap_stack(detail::overlap_stack_compare(end_field))
		{
		}
	
	public:
		variant_preprocessor(
			variant_preprocessor_delegate &delegate,
			libbio::vcf_reader &reader,
			vector_type const &reference,
			std::string const chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count
		):
			variant_preprocessor(
				delegate,
				reader,
				*variant_preprocessor::get_end_field_ptr(reader),
				reference,
				chr_name,
				donor_count,
				chr_count
			)
		{
		}
		
		variant_graph const &variant_graph() const { return m_graph; }
		
		void process(std::vector <std::string> const &field_names_for_filter_by_assigned);
		
	protected:
		void update_sample_names();
		void process_subgraph(std::size_t const prev_overlap_end);
		void calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx);
		void calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx, std::size_t const max_in_alt_edge_aligned_pos);
		void assign_alt_edge_labels_and_queue(libbio::variant const &var, std::size_t const node_idx, std::size_t const alt_edge_start_idx);
		void finalize_graph();
	};
}

#endif
