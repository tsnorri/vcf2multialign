/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/combine.hpp>
#include <vcf2multialign/alt_checker.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {
	
	bool alt_checker::check_alt_seq(std::string_view const &alt) const
	{
		for (auto const c : alt)
		{
			if (! ('A' == c || 'C' == c || 'G' == c || 'T' == c || 'N' == c))
				return false;
		}
		
		return true;
	}
	
	
	void alt_checker::check_variant(transient_variant const &var)
	{
		auto const lineno(var.lineno());
		auto &valid_alts(m_valid_alts_by_lineno.at(lineno));

		// Check that the alt sequence is something that can be handled.
		valid_alts.clear();
		valid_alts.reserve(var.alts().size());
		
		uint8_t i(0);
		
		if (sv_handling::DISCARD == m_sv_handling_method)
		{
			for (auto const &ref : boost::combine(var.alts(), var.alt_sv_types()))
			{
				++i;
				
				auto const alt_svt(ref.get <1>());
				if (sv_type::NONE != alt_svt)
					continue;
				
				auto const &alt(ref.get <0>());
				if (!check_alt_seq(alt))
				{
					m_error_logger->log_invalid_alt_seq(lineno, i, std::string(alt));
					continue;
				}
				
				valid_alts.emplace_back(i);
			}
		}
		else
		{
			for (auto const &ref : boost::combine(var.alts(), var.alt_sv_types()))
			{
				++i;
				
				auto const alt_svt(ref.get <1>());
				switch (alt_svt)
				{
					case sv_type::NONE:
					{
						auto const &alt(ref.get <0>());
						if (check_alt_seq(alt))
							valid_alts.emplace_back(i);
						else
							m_error_logger->log_invalid_alt_seq(lineno, i, std::string(alt));
						
						break;
					}
					
					case sv_type::DEL:
					case sv_type::DEL_ME:
						valid_alts.emplace_back(i);
						break;
					
					case sv_type::INS:
					case sv_type::DUP:
					case sv_type::INV:
					case sv_type::CNV:
					case sv_type::DUP_TANDEM:
					case sv_type::INS_ME:
					case sv_type::UNKNOWN:
						m_error_logger->log_skipped_structural_variant(lineno, i, alt_svt);
						break;
					
					default:
						fail("Unexpected structural variant type.");
						break;
				}
			}
		}
		
		std::sort(valid_alts.begin(), valid_alts.end());
	}
	
	
	void alt_checker::check_all_variants(vcf_reader &reader)
	{
		reader.reset();
		reader.set_parsed_fields(vcf_field::ALT);
		size_t last_position(0);
		bool should_continue(false);
		do {
			reader.fill_buffer();
			should_continue = reader.parse(
				[
					this,
					&last_position
				]
				(transient_variant const &var)
				-> bool
				{
					// Verify that the positions are in increasing order.
					auto const pos(var.zero_based_pos());
					
					always_assert(last_position <= pos, "Positions not in increasing order");
					check_variant(var);
					
					return true;
				}
			);
		} while (should_continue);
	}
}
