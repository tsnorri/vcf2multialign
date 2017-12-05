/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/util.hh>


namespace vcf2multialign {
	
	std::size_t strlen_utf8(std::string const &str)
	{
		std::size_t retval(0);
		for (auto const c : str)
		{
			if (0x80 != (0xc0 & c))
				++retval;
		}
		
		return retval;
	}
}
