/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <iomanip>
#include <sys/ioctl.h>
#include <vcf2multialign/progress_bar.hh>
#include <vcf2multialign/status_logger.hh>


namespace v2m = vcf2multialign;


namespace vcf2multialign {
	
	void status_logger::stop()
	{
		// Deallocating a suspended source is considered a bug.
		dispatch_resume(*m_message_timer);
		dispatch_source_cancel(*m_message_timer);
		dispatch_source_cancel(*m_signal_source);
	}
	
	
	void status_logger::finish_logging()
	{
		dispatch_suspend(*m_message_timer);
		std::cerr << std::endl;
	}
	
	
	void status_logger::log_message_counting(std::string const &message)
	{
		auto const message_len(strlen_utf8(message));
		m_message = message;
		m_message_length = message_len;

		v2m::dispatch(this).source_set_event_handler <&status_logger::update_count>(*m_message_timer);
		dispatch_resume(*m_message_timer);
	}
	
	
	void status_logger::log_message_progress_bar(std::string const &message)
	{
		auto const message_len(strlen_utf8(message));
		m_message = message;
		m_message_length = message_len;
		
		v2m::dispatch(this).source_set_event_handler <&status_logger::update_progress_bar>(*m_message_timer);
		dispatch_resume(*m_message_timer);
	}
	
	
	void status_logger::update_count()
	{
		auto fmt(
			boost::format("%s % d ")
			% m_message
			% boost::io::group(std::setw(m_window_width - m_message_length - 2), *m_record_count)
		);
		std::cerr << '\r' << fmt << std::flush;
	}
	
	
	void status_logger::update_progress_bar()
	{
		float const record_count(*m_record_count);
		float const current_record(*m_current_record);
		auto const half(m_window_width / 2);
		auto const pad(half < m_message_length ? 1 : half - m_message_length);
		auto const bar_width(m_window_width - half - 2);
		progress_bar(std::cerr, current_record / record_count, bar_width, pad, m_message);
	}


	void status_logger::handle_window_size_change()
	{
		struct winsize ws;
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws);
		m_window_width = ws.ws_col;
		if (!m_window_width)
			m_window_width = 80;
	}
}
