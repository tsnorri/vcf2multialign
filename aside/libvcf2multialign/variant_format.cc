/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/variant_format.hh>


namespace lb	= libbio;
namespace vcf	= libbio::vcf;


namespace vcf2multialign {
	
	void variant_format::reader_will_update_format(vcf::reader &reader)
	{
	}
	
	
	void variant_format::reader_did_update_format(vcf::reader &reader)
	{	
		this->assign_field_ptr("GT", gt);
	}
}
