/*
 * Copyright (c) 2019-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_VARIANT_GRAPH_GENERATOR_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_VARIANT_GRAPH_GENERATOR_HH

#include <deque>
#include <libbio/copyable_atomic.hh>
#include <libbio/int_matrix.hh>
#include <libbio/matrix.hh>
#include <libbio/vcf/subfield.hh>
#include <libbio/vcf/variant.hh>
#include <ostream>
#include <queue>
#include <string>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include <vcf2multialign/preprocess/path_sorted_variant.hh>
#include <vcf2multialign/preprocess/sample_indexer.hh>
#include <vcf2multialign/preprocess/sample_sorter.hh>
#include <vcf2multialign/preprocess/types.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_processor.hh>
#include <vcf2multialign/variant_processor_delegate.hh>
#include <vector>


namespace vcf2multialign {  namespace variant_graphs {
	class variant_graph_generator; // Fwd.
	
	struct variant_graph_generator_delegate
	{
		virtual ~variant_graph_generator_delegate() {}
		virtual void variant_graph_generator_will_handle_subgraph(libbio::vcf::variant const &first_var, std::size_t const variant_count, std::size_t const path_count) = 0;
	};
		
	struct variant_graph_single_pass_generator_delegate :	public virtual variant_graph_generator_delegate,
															public virtual variant_processor_delegate // For the case where processing was not done.
	{
	};
}}


namespace vcf2multialign { namespace variant_graphs { namespace detail {
	
	struct generator_node_description
	{
		std::size_t node_index{};
		std::size_t ref_position{};
		std::size_t alt_edge_start{};
		std::size_t alt_edge_count{};
	
		generator_node_description() = default;
	
		explicit generator_node_description(std::size_t ref_position_):
			ref_position(ref_position_)
		{
		}
	};
	
	
	// For having default implementations of the assignment operators.
	class variant_graph_generator_base
	{
	protected:
		typedef std::deque <libbio::vcf::variant>					variant_vector;
		typedef std::vector <detail::generator_node_description>	node_description_vector;
		typedef std::vector <std::size_t>							position_vector;
		
	protected:
		libbio::vcf::info_field_end const				*m_end_field{};
		variant_graph									m_graph;
		variant_vector									m_subgraph_variants;
		sample_indexer									m_sample_indexer;
		sample_sorter									m_sample_sorter;
		std::vector <std::string>						m_sample_names;
		node_description_vector							m_sorted_nodes;
		position_vector									m_start_positions;
		position_vector									m_end_positions;
		position_vector									m_end_positions_by_sample;
		std::vector <bool>								m_can_handle_alt;
		std::vector <std::uint16_t>						m_unhandled_alt_csum;
		std::size_t										m_output_lineno{};
		libbio::copyable_atomic <std::size_t>			m_processed_count{};
		
	public:
		variant_graph_generator_base() = default;
		
		variant_graph_generator_base(libbio::vcf::reader &reader, sample_sorter_delegate &delegate):
			m_end_field(reader.get_end_field_ptr()),
			m_sample_indexer(reader.current_variant()),
			m_sample_sorter(delegate, *m_end_field, m_sample_indexer)
		{
		}
	};
}}}


namespace vcf2multialign { namespace variant_graphs {
	
	// Generate a variant graph from a set of VCF records.
	class variant_graph_generator :	public detail::variant_graph_generator_base, sample_sorter_delegate
	{
	public:
		variant_graph_generator() = default;
		
		explicit variant_graph_generator(libbio::vcf::reader &reader):
			variant_graph_generator_base(reader, *this)
		{
		}
		
		class variant_graph const &variant_graph() const { return m_graph; }
		class variant_graph &variant_graph() { return m_graph; }
		class sample_sorter &sample_sorter() { return m_sample_sorter; }
		class sample_sorter const &sample_sorter() const { return m_sample_sorter; }
		
		inline void finish_copy_or_move(); // FIXME: should not be public.
		
		std::size_t processed_count() const { return m_processed_count.load(std::memory_order_relaxed); }
		
		virtual void sample_sorter_found_overlapping_variant(
			libbio::vcf::variant const &var,
			std::size_t const sample_idx,
			std::size_t const prev_end_pos
		) override {}
			
		virtual libbio::vcf::reader &vcf_reader() = 0;
		virtual vector_type const &reference() = 0;
		virtual variant_graph_generator_delegate &delegate() = 0;
		
	protected:
		void update_sample_names();
		void process_subgraph(std::size_t const prev_overlap_end);
		void combine_subgraph_variants_by_pos_and_ref();
		void reset_genotype_values_for_samples_with_overlapping_variants();
		std::tuple <std::size_t, std::size_t> create_subgraph_and_nodes();
		void generate_graph_setup();
		void finalize_graph();
	};
	
	
	class variant_graph_precalculated_generator final :	public variant_graph_generator
	{
	protected:
		variant_graph_generator_delegate				*m_delegate{};
		libbio::vcf::reader								*m_reader{};
		vector_type const								*m_reference{};
		preprocessing_result const						*m_preprocessing_result{};
		
	public:
		variant_graph_precalculated_generator() = default;
		
		variant_graph_precalculated_generator(
			variant_graph_generator_delegate &delegate,
			libbio::vcf::reader &reader,
			vector_type const &reference,
			preprocessing_result const &preprocessing_result
		):
			variant_graph_generator(reader),
			m_delegate(&delegate),
			m_reader(&reader),
			m_reference(&reference),
			m_preprocessing_result(&preprocessing_result)
		{
		}
		
		libbio::vcf::reader &vcf_reader() override { return *m_reader; }
		vector_type const &reference() override { return *m_reference; }
		variant_graph_generator_delegate &delegate() override { return *m_delegate; }
		
		void generate_graph(bool const should_start_from_current_variant);
	};
	
	
	class variant_graph_single_pass_generator final :	public variant_graph_generator,
														public variant_processor
	{
	protected:
		variant_graph_single_pass_generator_delegate	*m_delegate{};
		std::size_t										m_minimum_bridge_length{};
		
	public:
		variant_graph_single_pass_generator() = default;
		
		variant_graph_single_pass_generator(
			variant_graph_single_pass_generator_delegate &delegate,
			libbio::vcf::reader &reader,
			vector_type const &reference,
			std::string const &chr_name,
			std::size_t const minimum_bridge_length
		):
			variant_graph_generator(reader),
			variant_processor(reader, reference, chr_name),
			m_delegate(&delegate),
			m_minimum_bridge_length(minimum_bridge_length)
		{
		}
		
		libbio::vcf::reader &vcf_reader() override { return *this->m_reader; }
		vector_type const &reference() override { return *this->m_reference; }
		variant_graph_single_pass_generator_delegate &delegate() override { return *m_delegate; }
		
		void generate_graph(
			std::vector <std::string> const &field_names_for_filter_by_assigned,
			bool const should_start_from_current_variant
		);
	};
	
	
	void variant_graph_generator::finish_copy_or_move()
	{
		m_sample_sorter.set_delegate(*this);
		m_sample_sorter.set_sample_indexer(m_sample_indexer);
	}
}}

#endif
