/*
 * Copyright (c) 2017-2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/can_handle_variant_alts.hh>


namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {

	bool can_handle_variant_alts(vcf::transient_variant const &var)
	{
		// Keep.
		for (auto const &alt : var.alts())
		{
			auto const svt(alt.alt_sv_type);
			switch (svt)
			{
				// These structural variant types are currently handled.
				case vcf::sv_type::NONE:
				case vcf::sv_type::DEL:
				case vcf::sv_type::DEL_ME:
				case vcf::sv_type::UNKNOWN:
					return true;
					
				case vcf::sv_type::INS:
				case vcf::sv_type::DUP:
				case vcf::sv_type::INV:
				case vcf::sv_type::CNV:
				case vcf::sv_type::DUP_TANDEM:
				case vcf::sv_type::INS_ME:
				case vcf::sv_type::UNKNOWN_SV:
				default:
					break;
			}
		}
		
		return false;
	}
	
	
	bool can_handle_variant_alt(vcf::variant_alt_base const &alt)
	{
		auto const svt(alt.alt_sv_type);
		switch (svt)
		{
			case vcf::sv_type::NONE:
			case vcf::sv_type::DEL:
			case vcf::sv_type::DEL_ME:
			case vcf::sv_type::UNKNOWN:
				return true;
				
			case vcf::sv_type::INS:
			case vcf::sv_type::DUP:
			case vcf::sv_type::INV:
			case vcf::sv_type::CNV:
			case vcf::sv_type::DUP_TANDEM:
			case vcf::sv_type::INS_ME:
			case vcf::sv_type::UNKNOWN_SV:
			default:
				break;
		}

		return false;
	}
}
