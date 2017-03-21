/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTIL_HH
#define VCF2MULTIALIGN_UTIL_HH

#include <iostream>


namespace vcf2multialign {
	
	inline void fail(char const *message)
	{
		std::cerr << message << std::endl;
		abort();
	}
	
	
	inline void always_assert(bool check) { if (!check) abort(); }
	
	
	inline void always_assert(bool check, char const *message)
	{
		if (!check)
			fail(message);
	}
	
	
	template <typename t_fn>
	inline void always_assert(bool check, t_fn fn)
	{
		if (!check)
		{
			fn();
			abort();
		}
	}
}

#endif
