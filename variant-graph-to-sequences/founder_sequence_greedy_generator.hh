/*
 * Copyright (c) 2020–2022 Tuukka Norri
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
		typedef std::map <std::size_t, std::size_t>				removed_count_map;
		
	protected:
		std::size_t									m_founder_count{};
		std::size_t									m_tail_length{};
		bool										m_output_all_substrings{};		// Can be used to select the algorithm for choosing the substrings. For now, the greedy algorithm (that does not guarantee choosing all the substrings) is always enabled.
		bool										m_fill_unassigned_with_ref{};
		bool										m_should_remove_mid{};
		
	public:
		founder_sequence_greedy_generator(
			char const *subprocess_executable_path,
			std::size_t const founder_count,
			std::size_t const tail_length,
			bool const output_reference,
			bool const fill_unassigned_with_ref,
			bool const should_remove_mid,
			bool const may_overwrite
		):
			sequence_generator_base(subprocess_executable_path, output_reference, may_overwrite),
			m_founder_count(founder_count),
			m_tail_length(tail_length),
			m_fill_unassigned_with_ref(fill_unassigned_with_ref),
			m_should_remove_mid(should_remove_mid)
		{
		}
		
	protected:
		void do_work() override;
		
		void process_graph_and_output(
			output_stream_vector &output_files,
			removed_count_map &removed_counts,
			progress_indicator_delegate &progress_delegate
		) const;
		
		std::string output_path(std::size_t const file_idx) const override { return std::to_string(file_idx); }
	};
}

#endif
