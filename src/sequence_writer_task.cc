/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/tasks/sequence_writer_task.hh>


namespace vcf2multialign {
	
	void sequence_writer_task::execute()
	{
		m_vcf_reader.reset();
		m_vcf_reader.set_parsed_fields(vcf_field::ALT);
		std::size_t last_position(0);
		bool should_continue(false);
		do {
			m_vcf_reader.fill_buffer();
			should_continue = m_vcf_reader.parse(
				[
					this,
					&last_position
				]
				(transient_variant const &var)
				-> bool
				{
					// Verify that the positions are in increasing order.
					auto const pos(var.zero_based_pos());
					
					always_assert(last_position <= pos, "Positions not in increasing order");
					
					variant_stats::handle_variant(var);
					m_sequence_writer.handle_variant(var);
					
					return true;
				}
			);
		} while (should_continue);
		m_sequence_writer.finish();
		m_delegate->task_did_finish(*this);
	}
}
