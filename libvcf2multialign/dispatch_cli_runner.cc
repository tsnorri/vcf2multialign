/*
 * Copyright (c) 2019–2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/utility/dispatch_cli_runner.hh>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void dispatch_cli_runner::run(bool needs_progress_indicator)
	{
		try
		{
			if (needs_progress_indicator)
				install_progress_indicator();
			
			do_work();
			finish();
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			log_assertion_failure(exc);
		}
		catch (std::exception const &exc)
		{
			log_exception(exc);
		}
		catch (...)
		{
			log_unknown_exception();
		}
	}
	
	
	void dispatch_cli_runner::finish()
	{
		dispatch_async(dispatch_get_main_queue(), ^{
			m_progress_indicator.uninstall();
		});
	}
	
	
	void dispatch_cli_runner::log_assertion_failure(lb::assertion_failure_exception const &exc)
	{
		end_logging_no_update();
		lb::dispatch_async_fn(dispatch_get_main_queue(), [this, exc](){
			log_assertion_failure_exception(exc);
			m_progress_indicator.uninstall();
		});
	}
	
	
	void dispatch_cli_runner::log_exception(std::exception const &exc)
	{
		end_logging_no_update();
		std::string const reason(exc.what());
		lb::dispatch_async_fn(dispatch_get_main_queue(), [this, reason](){
			std::cerr << "Caught an exception: " << reason << '\n';
			m_progress_indicator.uninstall();
		});
	}
	
	
	void dispatch_cli_runner::log_unknown_exception()
	{
		end_logging_no_update();
		dispatch_async(dispatch_get_main_queue(), ^{
			std::cerr << "Caught an unknown exception.\n";
			m_progress_indicator.uninstall();
		});
	}
}
