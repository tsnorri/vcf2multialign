/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/logging_memory_resource.hh>


namespace vcf2multialign {
	
	void logging_memory_resource_base::log_status(status_logger &logger)
	{
		auto const bytes((std::size_t(m_allocated_bytes)));
		logger.log([this, bytes](){
			std::cerr << m_name << ": allocated bytes: " << bytes << std::endl;
		});
	}
}
