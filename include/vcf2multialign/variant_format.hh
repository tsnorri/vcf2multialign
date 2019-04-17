/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_FORMAT_HH
#define VCF2MULTIALIGN_VARIANT_FORMAT_HH

#include <libbio/vcf/variant.hh>
#include <libbio/vcf/variant_format.hh>
#include <libbio/vcf/vcf_subfield.hh>


namespace vcf2multialign {
	
	// Custom fields for 1000G variants.
	using vcf_info_field_cipos	= libbio::vcf_info_field <2, libbio::vcf_metadata_value_type::FLOAT>;
	using vcf_info_field_ciend	= libbio::vcf_info_field <2, libbio::vcf_metadata_value_type::FLOAT>;
	using vcf_info_field_svlen	= libbio::vcf_info_field <1, libbio::vcf_metadata_value_type::INTEGER>;
	using vcf_info_field_svtype	= libbio::vcf_info_field <1, libbio::vcf_metadata_value_type::STRING>;
	
	
	struct variant_format final : public libbio::variant_format
	{
		libbio::vcf_genotype_field_gt	*gt{};
		
		virtual void reader_will_update_format(libbio::vcf_reader &reader) override;
		virtual void reader_did_update_format(libbio::vcf_reader &reader) override;
	};
	
	
	inline variant_format const &get_variant_format(libbio::variant const &var)				{ return static_cast <variant_format const &>(var.get_format()); }
	inline variant_format const &get_variant_format(libbio::transient_variant const &var)	{ return static_cast <variant_format const &>(var.get_format()); }
}

#endif
