/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_BED_READER_DELEGATE_HH
#define VCF2MULTIALIGN_BED_READER_DELEGATE_HH

#include <libbio/vcf/region_variant_validator.hh>


namespace vcf2multialign {
	class bed_reader_delegate : public libbio::vcf::region_variant_validator_bed_reader_delegate
	{
		typedef libbio::vcf::region_variant_validator_bed_reader_delegate base_t;
		
	public:
		using base_t::base_t;
		
		void bed_reader_reported_error(std::size_t const lineno) override
		{
			std::cerr << "ERROR: Parse error in BED file on line " << lineno << ".\n";
			std::exit(EXIT_FAILURE);
		}
	};
}

#endif
