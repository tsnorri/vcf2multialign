/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_FOUNDER_SEQUENCE_GREEDY_GENERATOR_HH
#define VCF2MULTIALIGN_FOUNDER_SEQUENCE_GREEDY_GENERATOR_HH

#include "sequence_generator_base.hh"


namespace vcf2multialign {
	
	class progress_indicator_delegate; // Fwd
	
	
	class founder_sequence_greedy_generator final : public sequence_generator_base
	{
	public:
		typedef sequence_generator_base::output_stream_vector	output_stream_vector;
		typedef sequence_generator_base::stream_position_list	stream_position_list;
		
	protected:
		std::size_t									m_founder_count{};
		bool										m_replace_duplicates_with_n{};
		
	public:
		founder_sequence_greedy_generator(
			std::size_t const founder_count,
			bool const output_reference,
			bool const replace_duplicates_with_n,
			bool const may_overwrite
		):
			sequence_generator_base(output_reference, may_overwrite),
			m_founder_count(founder_count)
		{
		}
		
		void output_sequences() override;
		
	protected:
		void process_graph_and_output(output_stream_vector &output_files, progress_indicator_delegate &progress_delegate) const;
	};
}

#endif
