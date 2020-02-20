/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCFCONTAINS_VCF_RECORD_GENERATOR_HH
#define VCF2MULTIALIGN_VCFCONTAINS_VCF_RECORD_GENERATOR_HH

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/vcf_processor.hh>


namespace vcf2multialign {
	
	class vcf_record_generator : public vcf_processor
	{
	protected:
		libbio::vcf_reader::parser_state	m_parser_state;
		
	public:
		void prepare();
		bool next_variant(libbio::variant &out_var);
	};
}

#endif
