/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cmath>
#include <vcf2multialign/progress_bar.hh>

namespace v2m = vcf2multialign;


namespace vcf2multialign {
	void progress_bar(std::ostream &stream, float const value, std::size_t const length, std::size_t const pad, std::string const &title)
	{
		std::string const blocks[]{u8"▏", u8"▎", u8"▍", u8"▌", u8"▋", u8"▊", u8"▉", u8"█"};
		
		assert(0.0f <= value);
		assert(value <= 1.0f);
		
		auto const v(value * length);
		auto const integral(std::floor(v));
		auto const fraction(v - integral);
		
		stream << '\r' << title ;
		for (std::size_t i(0); i < pad; ++i)
			stream << ' ';
		
		stream << u8"▕";
		for (std::size_t i(0); i < integral; ++i)
			stream << u8"█";
		
		if (v != integral)
		{
			auto const block_idx(static_cast <std::size_t>(std::floor(fraction / 0.125)));
			stream << blocks[block_idx];

			for (std::size_t i(1 + integral); i < length; ++i)
				stream << ' ';
		}
		
		stream << u8"▏" << std::flush;
	}
}
