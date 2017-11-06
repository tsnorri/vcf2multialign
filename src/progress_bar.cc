/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <cmath>
#include <vcf2multialign/progress_bar.hh>


namespace {
	enum {
		PERCENT_WIDTH = 5,
		TIME_WIDTH = 9
	};
}


namespace ch = std::chrono;


namespace vcf2multialign {
	
	void progress_bar(
		std::ostream &stream,
		float const value,
		std::size_t length,
		std::size_t const pad,
		std::string const &title,
		ch::time_point <ch::steady_clock> start_time
	)
	{
		std::string const blocks[]{u8"▏", u8"▎", u8"▍", u8"▌", u8"▋", u8"▊", u8"▉", u8"█"};
		
		assert(0.0f <= value);
		assert(value <= 1.0f);
		
		if (PERCENT_WIDTH < length)
		{
			length -= PERCENT_WIDTH;
			
			stream << '\r' << title;
			for (std::size_t i(0); i < pad; ++i)
				stream << ' ';
			
			// Output the elapsed time.
			if (TIME_WIDTH < length)
			{
				length -= TIME_WIDTH;
				auto const current_time(ch::steady_clock::now());
				ch::duration <double> const elapsed_seconds(current_time - start_time);
				auto const seconds(ch::duration_cast <ch::seconds> (elapsed_seconds).count() % 60);
				auto const minutes(ch::duration_cast <ch::minutes> (elapsed_seconds).count() % 60);
				auto const hours(ch::duration_cast <ch::hours> (elapsed_seconds).count());
				stream << (boost::format("%02d:%02d:%02d ") % hours % minutes % seconds);
			}
			
			// Output the progress bar.
			auto const v(value * length);
			auto const integral(std::floor(v));
			auto const fraction(v - integral);
		
			//stream << u8"▕";
			for (std::size_t i(0); i < integral; ++i)
				stream << u8"█";
		
			if (v != integral)
			{
				// Divide the fraction by the block size, 1/8th.
				auto const block_idx(static_cast <std::size_t>(std::floor(fraction / 0.125)));
				stream << blocks[block_idx];

				for (std::size_t i(1 + integral); i < length; ++i)
					stream << ' ';
			}
		
			//stream << u8"▏";
		}
		
		// Display the value as percentage.
		stream << (boost::format("% 4.0f%%") % (100.0 * value)) << std::flush;
		
	}
}
