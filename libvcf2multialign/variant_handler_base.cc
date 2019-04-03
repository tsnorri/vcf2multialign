/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/combine.hpp>
#include <vcf2multialign/variant_handler_base.hh>


namespace lb = libbio;


namespace vcf2multialign {
	
	bool variant_handler_base::check_alt_seq(std::string const &alt) const
	{
		for (auto const c : alt)
		{
			if (! ('A' == c || 'C' == c || 'G' == c || 'T' == c || 'N' == c))
				return false;
		}
		
		return true;
	}
	
	
	// FIXME: check from the VCF spec if INS could be handled.
	void variant_handler_base::fill_valid_alts(lb::variant const &var)
	{
		// Check that the alt sequence is something that can be handled.
		m_valid_alts.clear();
		std::size_t i(0);
		auto const lineno(var.lineno());
		
		if (sv_handling::DISCARD == m_sv_handling_method)
		{
			for (auto const &ref : boost::combine(var.alts(), var.alt_sv_types()))
			{
				++i;
				
				auto const alt_svt(ref.get <1>());
				if (lb::sv_type::NONE != alt_svt)
					continue;
				
				auto const &alt(ref.get <0>());
				if (!check_alt_seq(alt))
				{
					m_error_logger->log_invalid_alt_seq(lineno, i, alt);
					continue;
				}
				
				m_valid_alts.emplace(i);
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
					case lb::sv_type::NONE:
					{
						auto const &alt(ref.get <0>());
						if (check_alt_seq(alt))
							m_valid_alts.emplace(i);
						else
							m_error_logger->log_invalid_alt_seq(lineno, i, alt);
						
						break;
					}
					
					case lb::sv_type::DEL:
					case lb::sv_type::DEL_ME:
						m_valid_alts.emplace(i);
						break;
					
					case lb::sv_type::INS:
					case lb::sv_type::DUP:
					case lb::sv_type::INV:
					case lb::sv_type::CNV:
					case lb::sv_type::DUP_TANDEM:
					case lb::sv_type::INS_ME:
					case lb::sv_type::UNKNOWN:
						m_error_logger->log_skipped_structural_variant(lineno, i, alt_svt);
						break;
					
					default:
						libbio_fail("Unexpected structural variant type.");
						break;
				}
			}
		}
	}
}
