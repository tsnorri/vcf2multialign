/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_PROGRESS_INDICATOR_MANAGER_HH
#define VCF2MULTIALIGN_UTILITY_PROGRESS_INDICATOR_MANAGER_HH

#include <libbio/progress_indicator.hh>


namespace vcf2multialign {

	// Maintain a progress indicator.
	class progress_indicator_manager
	{
	protected:
		libbio::progress_indicator			m_progress_indicator;

	public:
		libbio::progress_indicator &progress_indicator() { return m_progress_indicator; }

		void install_progress_indicator() { if (m_progress_indicator.is_stderr_interactive()) m_progress_indicator.install(); }
		void uninstall_progress_indicator() { m_progress_indicator.uninstall(); }

		void end_logging() { m_progress_indicator.end_logging(); }
		void end_logging_no_update() { m_progress_indicator.end_logging_no_update(); }
		void finish();
		void finish_mt();
		
		void log_assertion_failure_and_exit(libbio::assertion_failure_exception const &exc);
		void log_exception_and_exit(std::exception const &exc);
		void log_unknown_exception_and_exit();
	};
}

#endif
