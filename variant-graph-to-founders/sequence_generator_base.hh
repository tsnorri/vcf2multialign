/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_SEQUENCE_GENERATOR_BASE_HH
#define VCF2MULTIALIGN_SEQUENCE_GENERATOR_BASE_HH

#include <libbio/file_handling.hh>
#include <list>
#include <vcf2multialign/preprocess/variant_graph.hh>


namespace vcf2multialign {
	
	struct stream_position
	{
		std::size_t node{};
		std::size_t stream_number{};
		
		bool operator<(stream_position const &other) const { return node < other.node; }
	};
	
	
	class alt_edge_handler_base;
	
	
	class sequence_generator_base
	{
	public:
		typedef std::vector <libbio::file_ostream>	output_stream_vector;
		typedef std::list <stream_position>			stream_position_list;
		
	public:
		virtual ~sequence_generator_base() {}
		
		void output_sequences(
			char const *reference_path,
			char const *input_graph_path,
			char const *reference_seq_name,
			std::size_t const chunk_size,
			bool const output_reference,
			bool const may_overwrite
		);
		
	protected:
		virtual std::size_t const get_stream_count(variant_graph const &graph, bool const get_stream_count) const = 0;
		
		virtual std::unique_ptr <alt_edge_handler_base> make_alt_edge_handler(
			std::string_view const &reference_sv,
			variant_graph const &graph,
			output_stream_vector &output_files
		) const = 0;
		
		virtual void open_output_file(std::size_t const idx, libbio::file_ostream &of, libbio::writing_open_mode const mode, variant_graph const &graph) const = 0;
		
		void output_chunk(std::string_view const &reference_sv, variant_graph const &graph, output_stream_vector &output_files);
	};
	
	
	// Used a helper class b.c. I did not want to replace O(mn) inlineable matrix accesses with virtual function calls.
	// Not sure if it was worth the trouble, though.
	class alt_edge_handler_base
	{
	public:
		typedef sequence_generator_base::output_stream_vector	output_stream_vector;
		typedef sequence_generator_base::stream_position_list	stream_position_list;

	protected:
		std::string_view const					*m_reference_sv{};
		variant_graph const						*m_graph{};
		output_stream_vector					*m_output_files{};
		variant_graph::sample_path_vector const	*m_sample_paths{};
		variant_graph::path_edge_matrix const	*m_edges_by_path_and_variant{};
		std::size_t								m_lhs{};
		std::size_t								m_rhs{};
		std::size_t								m_alt_lhs{};
		std::size_t								m_ref_lhs{};
		std::size_t								m_ref_rhs{};
		std::size_t								m_aln_lhs{};
		std::size_t								m_aln_rhs{};

	public:
		alt_edge_handler_base() = default;
		
		alt_edge_handler_base(std::string_view const &reference_sv, variant_graph const &graph, output_stream_vector &output_files):
			m_reference_sv(&reference_sv),
			m_graph(&graph),
			m_output_files(&output_files)
		{
		}
		
		virtual ~alt_edge_handler_base() {}
		
		inline void set_subgraph_samples_and_edges(
			variant_graph::sample_path_vector const &sample_paths,
			variant_graph::path_edge_matrix const &edges_by_path_and_variant
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
		variant_graph::sample_path_vector const &sample_paths,
		variant_graph::path_edge_matrix const &edges_by_path_and_variant
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
