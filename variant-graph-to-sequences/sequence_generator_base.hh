/*
 * Copyright (c) 2019–2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_SEQUENCE_GENERATOR_BASE_HH
#define VCF2MULTIALIGN_SEQUENCE_GENERATOR_BASE_HH

#include <list>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include <vcf2multialign/utility/dispatch_cli_runner.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "sequence_output_handler.hh"


namespace vcf2multialign {
	
	struct stream_position
	{
		std::size_t node{};
		std::size_t stream_number{};
		
		bool operator<(stream_position const &other) const { return node < other.node; }
	};
	
	
	class alt_edge_handler_base;
	
	class sequence_generator_base :	public dispatch_cli_runner,
									public reference_controller
	{
	public:
		typedef vcf2multialign::output_stream_type	output_stream_type;
		typedef std::vector <output_stream_type>	output_stream_vector;
		typedef std::list <stream_position>			stream_position_list;
		
	protected:
		std::unique_ptr <sequence_output_handler>	m_output_handler{};
		variant_graphs::variant_graph				m_graph;
		bool										m_output_reference{};
		bool										m_may_overwrite{};
		
	public:
		sequence_generator_base(
			char const *subprocess_executable_path,
			bool const output_reference,
			bool const may_overwrite
		):
			m_output_handler(
				subprocess_executable_path
				? static_cast <sequence_output_handler *>(new subprocess_sequence_output_handler(subprocess_executable_path))
				: static_cast <sequence_output_handler *>(new file_sequence_output_handler)
			),
			m_output_reference(output_reference),
			m_may_overwrite(may_overwrite)
		{
		}
		
		virtual ~sequence_generator_base() {}
		
		void read_variant_graph(char const *variant_graph_path);
		virtual void prepare_variant_graph() {}
		variant_graphs::variant_graph const &variant_graph() const { return m_graph; }
		virtual std::string output_path(std::size_t const file_idx) const = 0;
	};
	
	
	class direct_matching_sequence_generator : public sequence_generator_base
	{
	public:
		typedef sequence_generator_base::output_stream_vector	output_stream_vector;
		typedef sequence_generator_base::stream_position_list	stream_position_list;
		
	protected:
		std::size_t									m_chunk_size{};
		
	public:
		direct_matching_sequence_generator(
			char const *subprocess_executable_path,
			std::size_t const chunk_size,
			bool const output_reference,
			bool const may_overwrite
		):
			sequence_generator_base(subprocess_executable_path, output_reference, may_overwrite),
			m_chunk_size(chunk_size)
		{
		}
		
	protected:
		void do_work() override;
		
		virtual std::size_t get_stream_count() const = 0;
		
		virtual std::unique_ptr <alt_edge_handler_base> make_alt_edge_handler(
			std::string_view const &reference_sv,
			variant_graphs::variant_graph const &graph,
			output_stream_vector &output_files
		) const = 0;
		
		void output_chunk(std::string_view const &reference_sv, output_stream_vector &output_files);
	};
	
	
	// Used a helper class b.c. I did not want to replace O(mn) inlineable matrix accesses with virtual function calls.
	// Not sure if it was worth the trouble, though.
	class alt_edge_handler_base
	{
	public:
		typedef sequence_generator_base::output_stream_vector	output_stream_vector;
		typedef sequence_generator_base::stream_position_list	stream_position_list;

	protected:
		std::string_view const									*m_reference_sv{};
		variant_graphs::variant_graph const						*m_graph{};
		output_stream_vector									*m_output_files{};
		variant_graphs::variant_graph::sample_path_vector const	*m_sample_paths{};
		variant_graphs::variant_graph::path_edge_matrix const	*m_edges_by_path_and_variant{};
		std::size_t												m_lhs{};
		std::size_t												m_rhs{};
		std::size_t												m_alt_lhs{};
		std::size_t												m_ref_lhs{};
		std::size_t												m_ref_rhs{};
		std::size_t												m_aln_lhs{};
		std::size_t												m_aln_rhs{};

	public:
		alt_edge_handler_base() = default;
		
		alt_edge_handler_base(
			std::string_view const &reference_sv,
			variant_graphs::variant_graph const &graph,
			output_stream_vector &output_files
		):
			m_reference_sv(&reference_sv),
			m_graph(&graph),
			m_output_files(&output_files)
		{
		}
		
		virtual ~alt_edge_handler_base() {}
		
		inline void set_subgraph_samples_and_edges(
			variant_graphs::variant_graph::sample_path_vector const &sample_paths,
			variant_graphs::variant_graph::path_edge_matrix const &edges_by_path_and_variant
		);
		
		inline void set_position(
			std::size_t const lhs,
			std::size_t const rhs,
			std::size_t const alt_lhs,
			std::size_t const ref_lhs,
			std::size_t const ref_rhs,
			std::size_t const aln_lhs,
			std::size_t const aln_rhs
		);
		
		virtual void handle_node_with_alt_edges(stream_position_list &files_waiting, std::size_t const subgraph_variant_idx) const = 0;
		
	protected:
		void handle_edge(stream_position &sp, std::size_t const edge_number) const;
	};
	
	
	void alt_edge_handler_base::set_subgraph_samples_and_edges(
		variant_graphs::variant_graph::sample_path_vector const &sample_paths,
		variant_graphs::variant_graph::path_edge_matrix const &edges_by_path_and_variant
	)
	{
		m_sample_paths = &sample_paths;
		m_edges_by_path_and_variant = &edges_by_path_and_variant;
	}
	
	
	void alt_edge_handler_base::set_position(
		std::size_t const lhs,
		std::size_t const rhs,
		std::size_t const alt_lhs,
		std::size_t const ref_lhs,
		std::size_t const ref_rhs,
		std::size_t const aln_lhs,
		std::size_t const aln_rhs
	)
	{
		m_lhs = lhs;
		m_rhs = rhs;
		m_alt_lhs = alt_lhs;
		m_ref_lhs = ref_lhs;
		m_ref_rhs = ref_rhs;
		m_aln_lhs = aln_lhs;
		m_aln_rhs = aln_rhs;
	}
}

#endif
