/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_LOGGING_MEMORY_RESOURCE_HH
#define VCF2MULTIALIGN_LOGGING_MEMORY_RESOURCE_HH

#include <atomic>
#include <vcf2multialign/status_logger.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	class logging_memory_resource_base
	{
	protected:
		std::atomic_size_t	m_allocated_bytes{};
		std::string			m_name{};
		
	public:
		logging_memory_resource_base(std::string const &name): m_name(name) {}
		virtual ~logging_memory_resource_base() {}
		void log_status(status_logger &logger);
	};
	
	
	template <typename t_resource>
	class logging_memory_resource : public logging_memory_resource_base
	{
	protected:
		typedef t_resource pool_resource_type;
		
	protected:
		pool_resource_type	m_pool_resource;
		std::atomic_size_t	m_allocated_bytes;
		
	public:
		using logging_memory_resource_base::logging_memory_resource_base;
		void *allocate(std::size_t bytes, std::size_t alignment) { m_allocated_bytes += bytes; return m_pool_resource.allocate(bytes, alignment); }
		void deallocate(void *p, std::size_t bytes, std::size_t alignment) { m_allocated_bytes -= bytes; m_pool_resource.deallocate(p, bytes, alignment); }
	};
}

#endif

