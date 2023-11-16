/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_SEQUENCE_WRITER_HH
#define VCF2MULTIALIGN_SEQUENCE_WRITER_HH

#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <vcf2multialign/variant_graph.hh>


namespace vcf2multialign {
	
	struct sequence_writing_delegate
	{
		typedef vcf2multialign::variant_graph	variant_graph;
		typedef variant_graph::sample_type		sample_type;
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::ploidy_type		ploidy_type;
		
		constexpr static inline auto const PLOIDY_MAX{variant_graph::PLOIDY_MAX};
		
		ploidy_type	chromosome_copy_index{PLOIDY_MAX};
		
		sequence_writing_delegate() = default;
		
		sequence_writing_delegate(ploidy_type const chromosome_copy_index_):
			chromosome_copy_index(chromosome_copy_index_)
		{
		}
		
		virtual ~sequence_writing_delegate() {}
		virtual void handle_node(variant_graph const &graph, node_type const node) = 0;
	};
	
	
	void output_sequence(
		sequence_type const &ref_seq,
		variant_graph const &graph,
		libbio::file_handle &fh,
		sequence_writing_delegate &delegate
	);
	
	
	void output_sequence(
		sequence_type const &ref_seq,
		variant_graph const &graph,
		libbio::file_ostream &stream,
		sequence_writing_delegate &delegate
	);
}

#endif
