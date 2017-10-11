/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PROGRESS_BAR_HH
#define VCF2MULTIALIGN_PROGRESS_BAR_HH

#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	void progress_bar(
		std::ostream &stream,
		float const value,
		std::size_t const length,
		std::size_t const pad,
		std::string const &title
	);
}

#endif
