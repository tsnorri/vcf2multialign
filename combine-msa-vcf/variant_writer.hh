/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VARIANT_WRITER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VARIANT_WRITER_HH

#include "output_handler.hh"


namespace vcf2multialign {
	
	class variant_writer : public output_handler
	{
	protected:
		std::string							m_output_chr_id;
		std::ostream						*m_os{};
		
	public:
		variant_writer() = default;
		
		variant_writer(std::ostream &os, std::string &&output_chr_id):
			m_output_chr_id(std::move(output_chr_id)),
			m_os(&os)
		{
		}
		
		void output_vcf_header() const;
		void handle_variant_description(variant_description &&desc) override;
	};
}

#endif
