/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/preprocess/variant_processor_delegate.hh>


namespace lb = libbio;


namespace vcf2multialign {
	
	void logging_variant_processor_delegate::variant_processor_no_field_for_identifier(std::string const &identifier)
	{
		std::cerr << "WARNING: Did not find a field for identifier “" << identifier << "”.\n";
	}
	
	
	void logging_variant_processor_delegate::variant_processor_found_variant_with_position_greater_than_reference_length(
		lb::transient_variant const &var
	)
	{
		std::cerr << "ERROR: Found a variant with a position greater than the reference length on line " << var.lineno() << "\n";
	}
	
	
	void logging_variant_processor_delegate::variant_processor_found_variant_with_no_suitable_alts(
		lb::transient_variant const &var
	)
	{
		std::cerr << "Line " << var.lineno() << ": Variant has no ALTs that could be handled.\n";
	}
		
	
	void logging_variant_processor_delegate::variant_processor_found_filtered_variant(
		lb::transient_variant const &var, lb::vcf_info_field_base const &field
	)
	{
		std::cerr << "Line " << var.lineno() << ": Variant has the field '" << field.get_metadata()->get_id() << "' set; skipping.\n";
	}
		
	
	void logging_variant_processor_delegate::variant_processor_found_variant_with_ref_mismatch(
		lb::transient_variant const &var, std::string_view const &ref_sub
	)
	{
		std::cerr << "WARNING: reference column mismatch on line " << var.lineno() << ": expected '" << ref_sub << "', got '" << var.ref() << "'\n";
	}
}