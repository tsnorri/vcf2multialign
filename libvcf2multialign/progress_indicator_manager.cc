/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/utility/progress_indicator_manager.hh>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void progress_indicator_manager::finish_mt()
	{
		m_progress_indicator.uninstall();
		std::exit(EXIT_SUCCESS);
	}
	
	
	void progress_indicator_manager::finish()
	{
		dispatch_async(dispatch_get_main_queue(), ^{
			finish_mt();
		});
	}
	
	
	void progress_indicator_manager::log_assertion_failure_and_exit(lb::assertion_failure_exception const &exc)
	{
		end_logging_no_update();
		lb::dispatch_async_fn(dispatch_get_main_queue(), [this, exc](){
			log_assertion_failure_exception(exc);
			m_progress_indicator.uninstall();
			std::exit(EXIT_FAILURE);
		});
	}
	
	
	void progress_indicator_manager::log_exception_and_exit(std::exception const &exc)
	{
		end_logging_no_update();
		std::string const reason(exc.what());
		lb::dispatch_async_fn(dispatch_get_main_queue(), [this, reason](){
			std::cerr << "Caught an exception: " << reason << '\n';
			m_progress_indicator.uninstall();
			std::exit(EXIT_FAILURE);
		});
	}
	
	
	void progress_indicator_manager::log_unknown_exception_and_exit()
	{
		end_logging_no_update();
		dispatch_async(dispatch_get_main_queue(), ^{
			std::cerr << "Caught an unknown exception.\n";
			m_progress_indicator.uninstall();
			std::exit(EXIT_FAILURE);
		});
	}
}
