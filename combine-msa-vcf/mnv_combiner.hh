/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_MNV_COMBINER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_MNV_COMBINER_HH

#include <optional>
#include "output_handler.hh"


namespace vcf2multialign {
	
	class mnv_combiner final : public output_handler
	{
	protected:
		std::optional <variant_description>	m_previous_desc;
		output_handler						*m_next_handler{};
	
	public:
		mnv_combiner() = default;
	
		mnv_combiner(output_handler &next_handler):
			m_next_handler(&next_handler)
		{
		}
		
		void handle_variant_description(variant_description &&desc) override;
		void finish();
	};
}

#endif
