/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/variant.hh>
#include "vcf_record_generator.hh"


namespace lb	= libbio;
namespace vcf	= libbio::vcf;


namespace vcf2multialign {

	void vcf_record_generator::prepare()
	{
		this->prepare_reader();
		
		auto &reader(this->m_vcf_reader);
		reader.reset();
	}
	
	
	bool vcf_record_generator::next_variant(vcf::variant &out_var)
	{
		bool retval(false);
		bool const st(this->m_vcf_reader.parse_one([this, &retval, &out_var](vcf::transient_variant const &var) {
			// Not reached on EOF.
			out_var = var;
			retval = true;
			return true;
		}, m_parser_state));
	
		return retval;
	}
}
