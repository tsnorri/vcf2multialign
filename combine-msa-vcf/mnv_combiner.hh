/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_MNV_COMBINER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_MNV_COMBINER_HH

#include <vector>
#include "output_handler.hh"


namespace vcf2multialign {
	
	class mnv_combiner final : public output_handler
	{
	protected:
		std::vector <variant_description>	m_previous_descs;
		output_handler						*m_next_handler{};
		std::uint16_t						m_ploidy{};
	
	public:
		mnv_combiner() = default;
	
		mnv_combiner(output_handler &next_handler, std::uint16_t const ploidy):
			m_next_handler(&next_handler),
			m_ploidy(ploidy)
		{
		}
		
		void handle_variant_description(variant_description &&desc) override;
		void finish();
	};
}

#endif
