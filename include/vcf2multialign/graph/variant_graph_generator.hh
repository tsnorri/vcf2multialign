/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GRAPH_VARIANT_PREPROCESSOR_HH
#define VCF2MULTIALIGN_GRAPH_VARIANT_PREPROCESSOR_HH

#include <libbio/copyable_atomic.hh>
#include <libbio/int_matrix.hh>
#include <libbio/matrix.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/vcf_subfield_def.hh>
#include <ostream>
#include <queue>
#include <string>
#include <vcf2multialign/graph/variant_graph.hh>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/sample_sorter.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_processor.hh>
#include <vcf2multialign/variant_processor_delegate.hh>
#include <vector>


namespace vcf2multialign {
	class variant_graph_generator; // Fwd.
	
	struct variant_graph_generator_delegate
	{
		virtual ~variant_graph_generator_delegate() {}
		virtual void variant_graph_generator_will_handle_subgraph(libbio::variant const &first_var, std::size_t const variant_count, std::size_t const path_count) = 0;
	};
		
	struct variant_graph_single_pass_generator_delegate :	public virtual variant_graph_generator_delegate,
															public virtual variant_processor_delegate // For the case where processing was not done.
	{
	};
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
		
		overlap_stack_compare() = default;
		
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
	
	class variant_graph_generator :	public sample_sorter_delegate
	{
	protected:
		typedef std::vector <libbio::variant>			variant_stack;
		typedef std::priority_queue <
			detail::overlap_stack_entry,
			std::vector <detail::overlap_stack_entry>,
			detail::overlap_stack_compare
		>												overlap_stack;
	
	protected:
		libbio::vcf_info_field_end const				*m_end_field{};
		variant_graph									m_graph;
		variant_stack									m_subgraph_variants;
		sample_indexer									m_sample_indexer;
		sample_sorter									m_sample_sorter;
		std::vector <std::string>						m_sample_names;
		overlap_stack									m_overlap_stack;
		std::size_t										m_minimum_subgraph_distance{};
		std::size_t										m_output_lineno{};
		libbio::copyable_atomic <std::size_t>			m_processed_count{};
	
	public:
		variant_graph_generator() = default;
		
		variant_graph_generator(
			libbio::vcf_reader &reader,
			std::size_t const donor_count,
			std::uint8_t const chr_count
		):
			m_end_field(reader.get_end_field_ptr()),
			m_sample_indexer(donor_count, chr_count), // donor_count and chr_count should be set in all cases.
			m_sample_sorter(*this, reader, m_sample_indexer),
			m_overlap_stack(detail::overlap_stack_compare(*m_end_field))
		{
		}
		
		class variant_graph const &variant_graph() const { return m_graph; }
		class variant_graph &variant_graph() { return m_graph; }
		class sample_sorter &sample_sorter() { return m_sample_sorter; }
		class sample_sorter const &sample_sorter() const { return m_sample_sorter; }
		class sample_indexer &sample_indexer() { return m_sample_indexer; }
		class sample_indexer const &sample_indexer() const { return m_sample_indexer; }
		
		std::size_t processed_count() const { return m_processed_count.load(std::memory_order_relaxed); }
		
		virtual void sample_sorter_found_overlapping_variant(
			libbio::variant const &var,
			std::size_t const sample_idx,
			std::size_t const prev_end_pos
		) override {}
			
		virtual libbio::vcf_reader &vcf_reader() = 0;
		virtual vector_type const &reference() = 0;
		virtual variant_graph_generator_delegate &delegate() = 0;
		
	protected:
		void update_sample_names();
		bool samples_have_alt(libbio::variant const &var, path_sorted_variant const &psv, std::size_t const given_alt_idx) const;
		void process_subgraph(std::size_t const prev_overlap_end);
		void calculate_aligned_ref_pos_for_new_node(std::size_t const node_idx);
		void update_aligned_ref_pos(std::size_t const node_idx, std::size_t const max_in_alt_edge_aligned_pos);
		void assign_alt_edge_labels_and_queue(libbio::variant const &var, std::size_t const node_idx, std::size_t const alt_edge_start_idx);
		void generate_graph_setup();
		void finalize_graph();
	};
	
	
	class variant_graph_precalculated_generator final :	public variant_graph_generator
	{
	protected:
		variant_graph_generator_delegate				*m_delegate{};
		libbio::vcf_reader								*m_reader{};
		vector_type const								*m_reference{};
		preprocessing_result const						*m_preprocessing_result{};
		
	public:
		variant_graph_precalculated_generator() = default;
		
		variant_graph_precalculated_generator(
			variant_graph_generator_delegate &delegate,
			libbio::vcf_reader &reader,
			vector_type const &reference,
			preprocessing_result const &preprocessing_result
		):
			variant_graph_generator(reader, preprocessing_result.donor_count, preprocessing_result.chr_count),
			m_delegate(&delegate),
			m_reader(&reader),
			m_reference(&reference),
			m_preprocessing_result(&preprocessing_result)
		{
		}
		
		libbio::vcf_reader &vcf_reader() override { return *m_reader; }
		vector_type const &reference() override { return *m_reference; }
		variant_graph_generator_delegate &delegate() override { return *m_delegate; }
		
		void generate_graph();
	};
	
	
	class variant_graph_single_pass_generator final :	public variant_graph_generator,
														public variant_processor
	{
	protected:
		variant_graph_single_pass_generator_delegate	*m_delegate{};
		
	public:
		variant_graph_single_pass_generator() = default;
		
		variant_graph_single_pass_generator(
			variant_graph_single_pass_generator_delegate &delegate,
			libbio::vcf_reader &reader,
			vector_type const &reference,
			std::string const &chr_name,
			std::size_t const donor_count,
			std::uint8_t const chr_count
		):
			variant_graph_generator(reader, donor_count, chr_count),
			variant_processor(reader, reference, chr_name),
			m_delegate(&delegate)
		{
		}
		
		libbio::vcf_reader &vcf_reader() override { return *this->m_reader; }
		vector_type const &reference() override { return *this->m_reference; }
		variant_graph_single_pass_generator_delegate &delegate() override { return *m_delegate; }
		
		void generate_graph(
			std::vector <std::string> const &field_names_for_filter_by_assigned
		);
	};
}

#endif
