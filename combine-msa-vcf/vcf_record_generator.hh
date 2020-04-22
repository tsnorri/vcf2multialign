/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VCF_RECORD_GENERATOR_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VCF_RECORD_GENERATOR_HH

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "types.hh"


namespace vcf2multialign {
	
	class vcf_record_generator : public vcf_processor
	{
	protected:
		libbio::vcf::reader::parser_state	m_parser_state;
		libbio::vcf::info_field_end const	*m_end_field{};
		
	public:
		void prepare();
		variant_record next_variant();
	};
}

#endif
