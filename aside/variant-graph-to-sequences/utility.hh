/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_TO_SEQUENCES_UTILITY_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_TO_SEQUENCES_UTILITY_HH

#include <atomic>
#include <libbio/dispatch.hh>


namespace vcf2multialign {
	
	class progress_indicator_delegate final : public libbio::progress_indicator_delegate
	{
	protected:
		std::size_t			m_progress_max{};
		std::atomic_size_t	m_progress_current{};
		
	public:
		progress_indicator_delegate(std::size_t const progress_max):
			m_progress_max(progress_max)
		{
		}
		
		virtual std::size_t progress_step_max() const { return m_progress_max; }
		virtual std::size_t progress_current_step() const { return m_progress_current.load(std::memory_order_relaxed); }
		virtual void progress_log_extra() const {}
		void advance() { m_progress_current.fetch_add(1, std::memory_order_relaxed); }
	};
	
	
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
}

#endif
