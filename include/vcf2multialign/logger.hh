/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_LOGGER_HH
#define VCF2MULTIALIGN_LOGGER_HH

#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/status_logger.hh>


namespace vcf2multialign {
	
	struct logger
	{
		class status_logger		status_logger;
		class error_logger		error_logger;
	};
}

#endif
