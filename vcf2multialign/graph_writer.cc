/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/haplotypes/graph_writer.hh>

namespace vcf2multialign
{
	// Dummy implementation.
	
	void graph_writer_impl::init(std::size_t n_rows)
	{
		graph.Init(n_rows);
	}
	
	void graph_writer_impl::process_segment(std::vector <std::vector <std::uint8_t>> const &segment_vec)
	{
		graph.AppendMatrixRange(segment_vec);
	}
	
	void graph_writer_impl::finish()
	{
		graph.CompleteConstruction();
		graph.PrintToGFA(*m_dst_file);
	}
}
