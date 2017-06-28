/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/types.hh>

namespace vcf2multialign {
	
	char const *to_string(sv_type const svt)
	{
		switch (svt)
		{
			case sv_type::NONE:
				return "(none)";
			
			case sv_type::DEL:
				return "DEL";
			
			case sv_type::INS:
				return "INS";
			
			case sv_type::DUP:
				return "DUP";
			
			case sv_type::INV:
				return "INV";
			
			case sv_type::CNV:
				return "CNV";
			
			case sv_type::DUP_TANDEM:
				return "DUP:TANDEM";
				
			case sv_type::DEL_ME:
				return "DEL:ME";
				
			case sv_type::INS_ME:
				return "INS:ME";
				
			case sv_type::UNKNOWN:
				return "(unknown)";
				
			default:
				return "(unexpected value)";
		}
	}
}
