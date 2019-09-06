/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_PROCESSOR_DELEGATE_HH
#define VCF2MULTIALIGN_VARIANT_PROCESSOR_DELEGATE_HH

#include <libbio/vcf/variant.hh>
#include <libbio/vcf/vcf_subfield_def.hh>
#include <string>
#include <string_view>


namespace vcf2multialign {
	
	struct variant_processor_delegate
	{
		virtual ~variant_processor_delegate() {}
		virtual void variant_processor_no_field_for_identifier(std::string const &identifier) = 0;
		virtual void variant_processor_found_variant_with_position_greater_than_reference_length(libbio::transient_variant const &var) = 0;
		virtual void variant_processor_found_variant_with_no_suitable_alts(libbio::transient_variant const &var) = 0;
		virtual void variant_processor_found_filtered_variant(libbio::transient_variant const &var, libbio::vcf_info_field_base const &field) = 0;
		virtual void variant_processor_found_variant_with_ref_mismatch(libbio::transient_variant const &var, std::string_view const &ref_sub) = 0;
	};
	
	
	struct logging_variant_processor_delegate : virtual public variant_processor_delegate
	{
		void variant_processor_no_field_for_identifier(std::string const &identifier) override;
		void variant_processor_found_variant_with_position_greater_than_reference_length(libbio::transient_variant const &var) override;
		void variant_processor_found_variant_with_no_suitable_alts(libbio::transient_variant const &var) override;
		void variant_processor_found_filtered_variant(libbio::transient_variant const &var, libbio::vcf_info_field_base const &field) override;
		void variant_processor_found_variant_with_ref_mismatch(libbio::transient_variant const &var, std::string_view const &ref_sub) override;
	};
}

#endif
