/*
 * Copyright (c) 2017-2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/can_handle_variant_alts.hh>


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {

	bool can_handle_variant_alts(lb::transient_variant const &var)
	{
		// Keep.
		for (auto const &alt : var.alts())
		{
			auto const svt(alt.alt_sv_type);
			switch (svt)
			{
				// These structural variant types are currently handled.
				case lb::sv_type::NONE:
				case lb::sv_type::DEL:
				case lb::sv_type::DEL_ME:
				case lb::sv_type::UNKNOWN:
					return true;
					
				case lb::sv_type::INS:
				case lb::sv_type::DUP:
				case lb::sv_type::INV:
				case lb::sv_type::CNV:
				case lb::sv_type::DUP_TANDEM:
				case lb::sv_type::INS_ME:
				case lb::sv_type::UNKNOWN_SV:
				default:
					break;
			}
		}
		
		return false;
	}
	
	
	bool can_handle_variant_alt(lb::variant_alt_base const &alt)
	{
		auto const svt(alt.alt_sv_type);
		switch (svt)
		{
			case lb::sv_type::NONE:
			case lb::sv_type::DEL:
			case lb::sv_type::DEL_ME:
			case lb::sv_type::UNKNOWN:
				return true;
				
			case lb::sv_type::INS:
			case lb::sv_type::DUP:
			case lb::sv_type::INV:
			case lb::sv_type::CNV:
			case lb::sv_type::DUP_TANDEM:
			case lb::sv_type::INS_ME:
			case lb::sv_type::UNKNOWN_SV:
			default:
				break;
		}

		return false;
	}
}
