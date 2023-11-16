/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/sequence_writer.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void output_sequence(
		sequence_type const &ref_seq,
		variant_graph const &graph,
		lb::file_ostream &stream,
		sequence_writing_delegate &delegate
	)
	{
		typedef variant_graph::position_type	position_type;
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::edge_type		edge_type;
		
		position_type ref_pos{};
		position_type aln_pos{};
		position_type next_ref_pos{};
		position_type next_aln_pos{};
		node_type current_node{};
		auto const limit(graph.node_count() - 1);
		while (current_node < limit)
		{
			delegate.handle_node(graph, current_node);
			
			std::size_t label_size{};
			if (sequence_writing_delegate::PLOIDY_MAX != delegate.chromosome_copy_index) // Always follow REF edges if outputting the aligned reference.
			{
				auto const &[edge_lb, edge_rb] = graph.edge_range_for_node(current_node);
				for (edge_type edge_idx(edge_lb); edge_idx < edge_rb; ++edge_idx)
				{
					if (graph.paths_by_chrom_copy_and_edge(edge_idx, delegate.chromosome_copy_index))
					{
						// Found an ALT edge to follow.
						auto const target_node(graph.alt_edge_targets[edge_idx]);
						auto const &label(graph.alt_edge_labels[edge_idx]);
						next_ref_pos = graph.reference_positions[target_node];
						next_aln_pos = graph.aligned_positions[target_node];
						libbio_assert_lte(label.size(), next_aln_pos - aln_pos);
						stream << label;
						current_node = target_node;
						label_size = label.size();
						goto continue_loop;
					}
				}
			}
			
			{
				next_ref_pos = graph.reference_positions[current_node + 1];
				next_aln_pos = graph.aligned_positions[current_node + 1];
				std::string_view const ref_part(ref_seq.data() + ref_pos, next_ref_pos - ref_pos);
				stream << ref_part;
				label_size = ref_part.size();
				++current_node;
			}
			
		continue_loop:
			std::fill_n(std::ostreambuf_iterator <char>(stream), next_aln_pos - aln_pos - label_size, '-');
			ref_pos = next_ref_pos;
			aln_pos = next_aln_pos;
		}
	}
	
	
	void output_sequence(
		sequence_type const &ref_seq,
		variant_graph const &graph,
		lb::file_handle &fh,
		sequence_writing_delegate &delegate
	)
	{
		lb::file_ostream stream;
		lb::open_stream_with_file_handle(stream, fh);
		output_sequence(ref_seq, graph, stream, delegate);
	}
}
