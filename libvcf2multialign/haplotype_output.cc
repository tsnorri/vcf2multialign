/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdint>
#include <libbio/file_handling.hh>
#include <ostream>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <sstream>
#include <vcf2multialign/output.hh>
#include <vcf2multialign/sequence_writer.hh>
#include <vcf2multialign/variant_graph.hh>

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

		if (m_should_output_reference)
		{
			std::stringstream fasta_id;
			if (m_chromosome_id)
				fasta_id << m_chromosome_id << '\t';
			fasta_id << "REF";

			::sequence_writing_delegate delegate;
			output_sequence(ref_seq, graph, stream, fasta_id.str().data(), m_should_output_unaligned, delegate);
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

				std::stringstream fasta_id;
				if (m_chromosome_id)
					fasta_id << m_chromosome_id << '\t';
				fasta_id << sample << '-' << (1 + chr_copy_idx);

				::sequence_writing_delegate delegate(graph, sample_idx, chr_copy_idx);
				output_sequence(ref_seq, graph, stream, fasta_id.str().data(), m_should_output_unaligned, delegate);
				stream << '\n';

				++seq_count;
				m_delegate->handled_sequences(seq_count);
			}
		}
	}


	void haplotype_output::output_separate(sequence_type const &ref_seq, variant_graph const &graph, bool const should_include_fasta_header)
	{
		typedef variant_graph::ploidy_type ploidy_type;

		if (m_should_output_reference)
		{
			// FIXME: Use std::format.
			std::stringstream dst_name;
			if (m_chromosome_id)
				dst_name << m_chromosome_id << '.';
			dst_name << "REF";
			if (should_include_fasta_header)
			{
				if (m_should_output_unaligned)
					dst_name << ".fa";
				else
					dst_name << ".a2m";
			}

			::sequence_writing_delegate delegate;
			output_sequence_file(ref_seq, graph, dst_name.str().data(), should_include_fasta_header, delegate);
		}

		for (auto const &[sample_idx, sample] : rsv::enumerate(graph.sample_names))
		{
			auto const ploidy(graph.sample_ploidy(sample_idx));
			for (auto const chr_copy_idx : rsv::iota(ploidy_type(0), ploidy))
			{
				m_delegate->will_handle_sample(sample, sample_idx, chr_copy_idx);

				// FIXME: Use std::format.
				std::stringstream dst_name;
				if (m_chromosome_id)
					dst_name << m_chromosome_id << '.';
				dst_name << sample << '.' << (1 + chr_copy_idx);
				if (should_include_fasta_header)
				{
					if (m_should_output_unaligned)
						dst_name << ".fa";
					else
						dst_name << ".a2m";
				}

				::sequence_writing_delegate delegate(graph, sample_idx, chr_copy_idx);
				output_sequence_file(ref_seq, graph, dst_name.str().data(), should_include_fasta_header, delegate);
			}
		}
	}
}
