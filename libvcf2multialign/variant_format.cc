/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/variant_format.hh>


namespace lb	= libbio;


#define ASSIGN_FIELD_PTR(ID, DST) do { \
	auto const it(this->m_fields_by_identifier.find(ID)); \
	libbio_always_assert_neq(it, m_fields_by_identifier.end()); \
	DST = dynamic_cast <decltype(DST)>(it->second.get()); \
} while (false)


namespace vcf2multialign {
	
	void variant_format::reader_will_update_format(lb::vcf_reader &reader)
	{
	}
	
	
	void variant_format::reader_did_update_format(lb::vcf_reader &reader)
	{	
		ASSIGN_FIELD_PTR("GT", gt);
	}
}
