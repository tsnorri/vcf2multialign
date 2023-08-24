/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VARIANT_FILTER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VARIANT_FILTER_HH

#include <ostream>
#include <string>
#include <vector>
#include "output_handler.hh"
#include "variant_description.hh"


namespace vcf2multialign {
	
	class variant_filter final : public output_handler
	{
	protected:
		std::vector <variant_description>	m_msa_output_variants;
		std::vector <variant_description>	m_vc_output_variants;
		output_handler						*m_next_handler{};
		std::size_t							m_current_pos{};
		
	protected:
		void clear_queue();
		
	public:
		variant_filter() = default;
		variant_filter(output_handler &next_handler):
			m_next_handler(&next_handler)
		{
		}
		
		void handle_variant_description(variant_description &&desc) override;
		void finish() override;
	};
}

#endif
