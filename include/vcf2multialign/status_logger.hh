/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_STATUS_LOGGER_HH
#define VCF2MULTIALIGN_STATUS_LOGGER_HH

#include <boost/format.hpp>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/util.hh>


namespace vcf2multialign {
	
	class status_logger
	{
	protected:
		dispatch_ptr <dispatch_source_t>	m_message_timer;
		dispatch_ptr <dispatch_source_t>	m_signal_source;
		std::string							m_message;
		std::atomic_size_t					*m_record_count{nullptr};
		std::atomic_size_t					*m_current_record{nullptr};
		std::size_t							m_window_width{0};
		std::size_t							m_message_length{0};
		
	public:
		status_logger(
			dispatch_ptr <dispatch_queue_t> &logging_queue,
			std::atomic_size_t *record_count,
			std::atomic_size_t *current_record
		):
			m_message_timer(
				dispatch_source_create(DISPATCH_SOURCE_TYPE_TIMER, 0, 0, *logging_queue),
				false
			),
			m_signal_source(
				dispatch_source_create(DISPATCH_SOURCE_TYPE_SIGNAL, SIGWINCH, 0, *logging_queue),
				false
			),
			m_record_count(record_count),
			m_current_record(current_record)
		{
			dispatch_source_set_timer(*m_message_timer, dispatch_time(DISPATCH_TIME_NOW, 0), 100000000, 10000000);
			dispatch(this).source_set_event_handler <&status_logger::handle_window_size_change>(*m_signal_source);
			dispatch(this).async <&status_logger::handle_window_size_change>(*logging_queue);
			
			dispatch_resume(*m_signal_source);
		}
		
		void stop();
		void finish_logging();
		void log_message_counting(std::string const &message);
		void log_message_progress_bar(std::string const &message);
		void update_count();
		void update_progress_bar();
		
	protected:
		void handle_window_size_change();
	};
}

#endif
