/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GRAPH_WRITER_HH
#define VCF2MULTIALIGN_GRAPH_WRITER_HH

#include <cstdint>
#include <iostream>
#include <vector>
#include <msa2dag/graph.h>


namespace vcf2multialign {
	
	class graph_writer_impl {
	protected:
		std::ostream		*m_dst_file{nullptr};
	
	public:
		graph_writer_impl() = default;
		graph_writer_impl(std::ostream &os): m_dst_file(&os) {}
	
		// Called before processing any segment.
		void init(std::size_t n_rows);
	
		// Process a segment that consists of std::vector <std::uint8_t>s of equal length by writing them
		// to m_dst_file as a graph.
		void process_segment(std::vector <std::vector <std::uint8_t>> const &segment_vec);
	
		// Called when process_segment will not be called again.
		void finish();
	private:
		MSA2DAG::Graph graph;
	};


	template <typename t_delegate>
	class graph_writer : public graph_writer_impl
	{
	protected:
		t_delegate			*m_delegate{nullptr};
	
	public:
		graph_writer() = default;
		graph_writer(t_delegate &delegate, std::ostream &os):
			graph_writer_impl(os),
			m_delegate(&delegate)
		{
		}
	
		void process_segment(std::vector <std::vector <std::uint8_t>> const &segment_vec)
		{
			graph_writer_impl::process_segment(segment_vec);
			m_delegate->graph_writer_did_process_segment(*this);
		}
	
		void finish()
		{
			graph_writer_impl::finish();
		}
	};
}

#endif
