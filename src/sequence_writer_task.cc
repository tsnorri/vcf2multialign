/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/tasks/sequence_writer_task.hh>


namespace vcf2multialign {
	
	void sequence_writer_task::prepare(class vcf_reader &reader)
	{
		reader.set_parsed_fields(vcf_field::ALT);
	}
	
	
	void sequence_writer_task::execute()
	{
		m_variant_handler.process_variants();
	}
	
	
	void sequence_writer_task::finish()
	{
		m_sequence_writer.finish();
		m_delegate->task_did_finish(*this);
	}
	
	
	void sequence_writer_task::handle_variant(variant &var)
	{
		variant_stats::handle_variant(var);
		m_sequence_writer.handle_variant(var);
	}
}
