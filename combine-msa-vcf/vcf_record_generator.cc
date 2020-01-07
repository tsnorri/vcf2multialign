/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/variant.hh>
#include "vcf_record_generator.hh"


namespace lb = libbio;


namespace vcf2multialign {

	void vcf_record_generator::prepare()
	{
		this->prepare_reader();
		
		auto &reader(this->m_vcf_reader);
		
		// Init the end field description.
		m_end_field = reader.get_end_field_ptr();
		
		reader.reset();
	}
	
	
	variant_record vcf_record_generator::next_variant()
	{
		variant_record retval;
		bool const st(this->m_vcf_reader.parse_one([this, &retval](lb::transient_variant const &var) {
			// Not reached on EOF.
			libbio_assert(m_end_field);
			retval.variant = var;
			retval.aligned_position = var.pos();
			retval.size = lb::variant_end_pos(var, *m_end_field) - var.zero_based_pos(); // FIXME: what do we store?
			return true;
		}, m_parser_state));
	
		return retval;
	}
}
