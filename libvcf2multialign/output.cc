/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/file_handling.hh>
#include <vcf2multialign/output.hh>
#include <vcf2multialign/sequence_writer.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void output::output_sequence_file(
		sequence_type const &ref_seq,
		variant_graph const &graph,
		char const * const dst_name,
		sequence_writing_delegate &delegate
	)
	{
		if (m_pipe_cmd)
		{
			auto proc(subprocess_type::subprocess_with_arguments({m_pipe_cmd, dst_name}));
			auto &fh(proc.stdin_handle());
			output_sequence(ref_seq, graph, fh, delegate);
			m_delegate->exit_subprocess(proc);
		}
		else
		{
			lb::file_handle fh(lb::open_file_for_writing(dst_name, lb::writing_open_mode::CREATE));
			output_sequence(ref_seq, graph, fh, delegate);
		}
	}
	
	
	void output::output_a2m(sequence_type const &ref_seq, variant_graph const &graph, char const * const dst_name)
	{
		if (m_pipe_cmd)
		{
			auto proc(subprocess_type::subprocess_with_arguments({m_pipe_cmd, dst_name}));
			auto &fh(proc.stdin_handle());
			
			{
				lb::file_ostream stream;
				lb::open_stream_with_file_handle(stream, fh);
				output_a2m(ref_seq, graph, stream);
			}
			
			m_delegate->exit_subprocess(proc);
		}
		else
		{
			lb::file_handle fh(lb::open_file_for_writing(dst_name, lb::writing_open_mode::CREATE));
			lb::file_ostream stream;
			lb::open_stream_with_file_handle(stream, fh);
			output_a2m(ref_seq, graph, stream);
		}
	}
}
