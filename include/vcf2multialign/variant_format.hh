/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_FORMAT_HH
#define VCF2MULTIALIGN_VARIANT_FORMAT_HH

#include <libbio/vcf/subfield.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/variant_format.hh>


namespace vcf2multialign {
	
	// Custom fields for 1000G variants.
	using vcf_info_field_cipos	= libbio::vcf::info_field <libbio::vcf::metadata_value_type::FLOAT,		2>;
	using vcf_info_field_ciend	= libbio::vcf::info_field <libbio::vcf::metadata_value_type::FLOAT,		2>;
	using vcf_info_field_svlen	= libbio::vcf::info_field <libbio::vcf::metadata_value_type::INTEGER,	1>;
	using vcf_info_field_svtype	= libbio::vcf::info_field <libbio::vcf::metadata_value_type::STRING,	1>;
	
	
	struct variant_format final : public libbio::vcf::variant_format
	{
		libbio::vcf::genotype_field_gt	*gt{};
		
		// Return a new empty instance of this class.
		virtual variant_format *new_instance() const override { return new variant_format(); }
		
		virtual void reader_will_update_format(libbio::vcf::reader &reader) override;
		virtual void reader_did_update_format(libbio::vcf::reader &reader) override;
	};
	
	
	inline variant_format const &get_variant_format(libbio::vcf::variant const &var)
	{
		libbio_assert(var.reader()->has_assigned_variant_format());
		return static_cast <variant_format const &>(var.get_format());
	}
	
	
	inline variant_format const &get_variant_format(libbio::vcf::transient_variant const &var)
	{
		libbio_assert(var.reader()->has_assigned_variant_format());
		return static_cast <variant_format const &>(var.get_format());
	}
}

#endif
