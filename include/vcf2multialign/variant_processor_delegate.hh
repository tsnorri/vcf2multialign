/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_PROCESSSOR_DELEGATE_HH
#define VCF2MULTIALIGN_VARIANT_PROCESSSOR_DELEGATE_HH

#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class variant;
	
	
	struct variant_processor_delegate
	{
		virtual ~variant_processor_delegate() {}
		
		virtual bool is_valid_alt(std::size_t const alt_idx) const = 0;
		
		virtual void enumerate_genotype(
			variant &var,
			std::size_t const sample_no,
			std::function <void(uint8_t, std::size_t, bool)> const &cb
		) = 0;
			
		virtual void assigned_alt_to_sequence(std::size_t const alt_idx) = 0;
		virtual void found_overlapping_alt(
			std::size_t const lineno,
			uint8_t const alt_idx,
			std::size_t const sample_no,
			uint8_t const chr_idx
		) = 0;
		virtual void handled_alt(std::size_t const alt_idx) = 0;
		virtual void handled_haplotypes(variant &var) = 0;
	};
}

#endif
