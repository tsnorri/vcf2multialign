/*
 * Copyright (c) 2019–2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_DISPATCH_CLI_RUNNER_HH
#define VCF2MULTIALIGN_UTILITY_DISPATCH_CLI_RUNNER_HH

#include <libbio/progress_indicator.hh>


namespace vcf2multialign {

	// Helper class for running a command line tool and maintaining a progress indicator.
	class dispatch_cli_runner
	{
	protected:
		libbio::progress_indicator			m_progress_indicator;

	public:
		libbio::progress_indicator &progress_indicator() { return m_progress_indicator; }
		
		void run(bool needs_progress_indicator);
		virtual void do_work() = 0;
		
		void uninstall_progress_indicator() { m_progress_indicator.uninstall(); }
		
		void end_logging() { m_progress_indicator.end_logging(); }
		void end_logging_no_update() { m_progress_indicator.end_logging_no_update(); }
		
	protected:
		void install_progress_indicator() { if (m_progress_indicator.is_stderr_interactive()) m_progress_indicator.install(); }
		void finish();
		void finish_mt() { m_progress_indicator.uninstall(); }
		
		void log_assertion_failure(libbio::assertion_failure_exception const &exc);
		void log_exception(std::exception const &exc);
		void log_unknown_exception();
	};
}

#endif
