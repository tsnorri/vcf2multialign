/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CXX14COMPAT_HH
#define VCF2MULTIALIGN_CXX14COMPAT_HH

#include <experimental/string_view>


#if __cplusplus < 201402L
namespace std {
	template <typename t_input_it_1, typename t_input_it_2>
	bool equal(t_input_it_1 it_1, t_input_it_1 const end_1, t_input_it_2 it_2, t_input_it_2 const end_2)
	{
		while (it_1 != end_1 && it_2 != end_2)
		{
			if (! (*it_1 == *it_2))
				return false;

			++it_1;
			++it_2;
		}

		// Ranges are equal iff they are of equal length and all elements compared equal.
		return (it_1 == end_1 && it_2 == end_2);
	}
}
#endif


#if __cplusplus < 201703L
// XXX Hack.
namespace std {
	using std::experimental::string_view;
	using std::experimental::optional;
	using std::experimental::nullopt;
}
#endif


#endif
