/*
 * Copyright (c) 2019-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/parse_error.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/variant_end_pos.hh>
#include "vcf_record_generator.hh"


namespace lb	= libbio;
namespace vcf	= libbio::vcf;


namespace vcf2multialign {

	void vcf_record_generator::prepare()
	{
		this->prepare_reader();
		
		auto &reader(this->m_vcf_reader);
		
		// Init the end field description.
		m_end_field = reader.get_end_field_ptr();
	}
	
	
	variant_record vcf_record_generator::next_variant()
	{
		variant_record retval;
		bool should_continue{true};
		bool have_input{true};
		
		try
		{
			while (should_continue && have_input)
			{
				have_input = this->m_vcf_reader.parse_one(
					[
						this,
						&should_continue,
						&retval
					](vcf::transient_variant const &var) {
						// Not reached on EOF.
					
						if (var.chrom_id() != m_chr_id)
							return true;
					
						should_continue = false; // Found a suitable variant.
						libbio_assert(m_end_field);
						retval.variant = var;
						retval.aligned_position = var.pos();
						retval.size = vcf::variant_end_pos(var, *m_end_field) - var.zero_based_pos(); // FIXME: what do we store?
						return true;
					},
					m_parser_state
				);
			}
		}
		catch (vcf::parse_error const &exc)
		{
			std::cerr << "ERROR: " << exc.what() << '\n';
			auto const value(exc.value());
			auto const field(exc.get_metadata());
			if (value)
				std::cerr << "Value: “" << (*value) << "”\n";
			if (field)
				std::cerr << "Field: " << (*field) << '\n';
			std::cerr << "Last characters from the buffer: " << this->m_vcf_reader.buffer_tail() << '\n';
			std::exit(EXIT_FAILURE);
		}
	
		return retval;
	}
}
