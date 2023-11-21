/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/file_handling.hh>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <type_traits>
#include <vcf2multialign/output.hh>
#include <vcf2multialign/sequence_writer.hh>

namespace ios	= boost::iostreams;
namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace v2m	= vcf2multialign;


namespace {
	
	struct sequence_writing_delegate final : public v2m::sequence_writing_delegate
	{
		virtual void handle_node(variant_graph const &graph, node_type const node) override {}
		
		using v2m::sequence_writing_delegate::sequence_writing_delegate;
		
		sequence_writing_delegate(variant_graph const &graph, sample_type const sample_idx, ploidy_type const chr_copy_idx):
			v2m::sequence_writing_delegate(graph.ploidy_csum[sample_idx] + chr_copy_idx)
		{
		}
	};
}


namespace vcf2multialign {
	
	void haplotype_output::output_a2m(
		sequence_type const &ref_seq,
		variant_graph const &graph,
		std::ostream &stream
	)
	{
		typedef variant_graph::ploidy_type	ploidy_type;
		
		std::uint32_t seq_count{1};
		
		{
			stream << ">REF\n";
			::sequence_writing_delegate delegate;
			output_sequence(ref_seq, graph, stream, delegate);
			stream << '\n';
			m_delegate->handled_sequences(seq_count);
		}
		
		auto const total_seq_count(graph.total_chromosome_copies());
		for (auto const &[sample_idx, sample] : rsv::enumerate(graph.sample_names))
		{
			auto const ploidy(graph.sample_ploidy(sample_idx));
			for (auto const chr_copy_idx : rsv::iota(ploidy_type(0), ploidy))
			{
				m_delegate->will_handle_sample(sample, sample_idx, chr_copy_idx);
				
				stream << '>' << sample << '-' << chr_copy_idx << '\n';
				::sequence_writing_delegate delegate(graph, sample_idx, chr_copy_idx);
				output_sequence(ref_seq, graph, stream, delegate);
				stream << '\n';
				
				++seq_count;
				m_delegate->handled_sequences(seq_count);
			}
		}
	}
	
	
	void haplotype_output::output_separate(sequence_type const &ref_seq, variant_graph const &graph)
	{
		typedef variant_graph::ploidy_type ploidy_type;
		
		{
			::sequence_writing_delegate delegate;
			output_sequence_file(ref_seq, graph, "REF", delegate);
		}
		
		for (auto const &[sample_idx, sample] : rsv::enumerate(graph.sample_names))
		{
			auto const ploidy(graph.sample_ploidy(sample_idx));
			for (auto const chr_copy_idx : rsv::iota(ploidy_type(0), ploidy))
			{
				m_delegate->will_handle_sample(sample, sample_idx, chr_copy_idx);

				// FIXME: Use std::format.
				std::stringstream dst_name;
				dst_name << sample;
				dst_name << '-';
				dst_name << chr_copy_idx;
				
				::sequence_writing_delegate delegate(graph, sample_idx, chr_copy_idx);
				output_sequence_file(ref_seq, graph, dst_name.str().data(), delegate);
			}
		}
	}
}
