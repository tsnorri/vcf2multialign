/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/can_handle_variant_alts.hh>
#include <vcf2multialign/variant_processor.hh>


namespace lb = libbio;


namespace vcf2multialign {
	
	void variant_processor::fill_filter_by_assigned(
		std::vector <std::string> const &field_names_for_filter_by_assigned,
		vcf_info_field_vector &filter_by_assigned,
		variant_processor_delegate &delegate
	)
	{
		auto const &fields(m_reader->info_fields());
		for (auto const &name : field_names_for_filter_by_assigned)
		{
			auto const it(fields.find(name));
			if (fields.end() == it)
			{
				delegate.variant_processor_no_field_for_identifier(name);
				continue;
			}
			
			filter_by_assigned.emplace_back(it->second.get());
		}
	}
	
	
	auto variant_processor::check_variant(
		lb::transient_variant const &var,
		vcf_info_field_vector const &filter_by_assigned,
		variant_processor_delegate &delegate
	) -> variant_check_status
	{
		// Check the chromosome name.
		if (var.chrom_id() != *m_chromosome_name)
		{
			delegate.variant_processor_found_variant_with_chrom_id_mismatch(var);
			return variant_check_status::ERROR;
		}
		
		auto const var_pos(var.zero_based_pos());
		if (! (var_pos < m_reference->size()))
		{
			delegate.variant_processor_found_variant_with_position_greater_than_reference_length(var);
			return variant_check_status::FATAL_ERROR;
		}
		
		if (!can_handle_variant_alts(var))
		{
			delegate.variant_processor_found_variant_with_no_suitable_alts(var);
			return variant_check_status::ERROR;
		}
		
		// Filter.
		for (auto const *field_ptr : filter_by_assigned)
		{
			if (field_ptr->has_value(var))
			{
				delegate.variant_processor_found_filtered_variant(var, *field_ptr);
				return variant_check_status::ERROR;
			}
		}
		
		// Compare the REF column against the reference sequence.
		{
			auto const &ref_col(var.ref());
			std::string_view const ref_sub(m_reference->data() + var_pos, ref_col.size());
			if (ref_col != ref_sub)
			{
				delegate.variant_processor_found_variant_with_ref_mismatch(var, ref_sub);
				return variant_check_status::ERROR;
			}
		}
		
		return variant_check_status::PASS;
	}
}
