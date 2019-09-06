/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_LOG_ASSERTION_FAILURE_HH
#define VCF2MULTIALIGN_UTILITY_LOG_ASSERTION_FAILURE_HH

#include <libbio/assert.hh>


namespace vcf2multialign {
	void log_assertion_failure_exception(libbio::assertion_failure_exception const &exc);
}

#endif
