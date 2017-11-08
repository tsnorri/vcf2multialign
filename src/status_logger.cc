/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <iomanip>
#include <sys/ioctl.h>
#include <vcf2multialign/progress_bar.hh>
#include <vcf2multialign/status_logger.hh>


namespace ch = std::chrono;
namespace v2m = vcf2multialign;


namespace vcf2multialign {
	
	void status_logger::install()
	{
		dispatch_queue_t main_queue(dispatch_get_main_queue());
		
		{
			dispatch_ptr <dispatch_source_t> message_timer(
				dispatch_source_create(DISPATCH_SOURCE_TYPE_TIMER, 0, 0, main_queue),
				false
			);
			
			dispatch_ptr <dispatch_source_t> signal_source(
				dispatch_source_create(DISPATCH_SOURCE_TYPE_SIGNAL, SIGWINCH, 0, main_queue),
				false
			);
		
			m_message_timer = std::move(message_timer);
			m_signal_source = std::move(signal_source);
		}
		
		dispatch_source_set_timer(*m_message_timer, dispatch_time(DISPATCH_TIME_NOW, 0), 333 * 1000 * 1000, 33 * 1000 * 1000);
		dispatch(this).source_set_event_handler <&status_logger::handle_window_size_change_mt>(*m_signal_source);
		dispatch(this).async <&status_logger::handle_window_size_change_mt>(main_queue);
		
		dispatch_resume(*m_signal_source);
	}
	
	
	void status_logger::uninstall()
	{
		// Deallocating a suspended source is considered a bug.
		dispatch_resume(*m_message_timer);
		
		dispatch_source_cancel(*m_message_timer);
		dispatch_source_cancel(*m_signal_source);
	}
	
	
	void status_logger::finish_logging()
	{
		dispatch(this).sync <&status_logger::finish_logging_mt>(dispatch_get_main_queue());
	}
	
	
	void status_logger::finish_logging_mt()
	{
		update_mt();
		dispatch_suspend(*m_message_timer);
		m_indicator_type = none;
		m_need_clear_line = false;
		m_timer_active = false;
		std::cerr << std::endl;
	}
	
	
	void status_logger::set_message(std::string const &message)
	{
		auto const message_len(strlen_utf8(message));
		std::lock_guard <std::mutex> guard(m_message_mutex);
		m_message = message;
		m_message_length = message_len;
	}
	
	
	void status_logger::log_message_counting(std::string const &message)
	{
		auto const message_len(strlen_utf8(message));
		
		{
			std::lock_guard <std::mutex> guard(m_message_mutex);
			m_message = message;
			m_message_length = message_len;
			m_indicator_type = counter;
			m_start_time = ch::steady_clock::now();
		}

		dispatch(this).async <&status_logger::resume_timer>(dispatch_get_main_queue());
	}
	
	
	void status_logger::log_message_progress_bar(std::string const &message)
	{
		auto const message_len(strlen_utf8(message));
		
		{
			std::lock_guard <std::mutex> guard(m_message_mutex);
			m_message = message;
			m_message_length = message_len;
			m_indicator_type = progress_bar;
			m_start_time = ch::steady_clock::now();
		}
		
		dispatch(this).async <&status_logger::resume_timer>(dispatch_get_main_queue());
	}
	
	
	void status_logger::update()
	{
		dispatch(this).async <&status_logger::update_mt>(dispatch_get_main_queue());
	}
	
	
	void status_logger::update_mt()
	{
		if (!m_timer_active)
			return;
		
		switch (m_indicator_type)
		{
			case counter:
			{
				std::lock_guard <std::mutex> guard(m_message_mutex);
				m_need_clear_line = true;
				auto const current_step(m_delegate->current_step());
				
				auto fmt(
					boost::format("%s % d ")
					% m_message
					% boost::io::group(std::setw(m_window_width - m_message_length - 2), current_step)
				);
				//std::cerr << '\r' << fmt << std::flush;
				std::cerr << "\33[2K\r" << fmt << std::flush;
				break;
			}
			
			case progress_bar:
			{
				std::lock_guard <std::mutex> guard(m_message_mutex);
				m_need_clear_line = true;
				
				float const step_count(m_delegate->step_count());
				float const current_step(m_delegate->current_step());
				float const progress_val(current_step / step_count);
				auto const half(m_window_width / 2);
				auto const pad(half < m_message_length ? 1 : half - m_message_length);

				v2m::progress_bar(std::cerr, progress_val, m_window_width - half - 1, pad, m_message, m_start_time);
			}
			
			case none:
			default:
				break;
		}
	}
	
	
	void status_logger::resume_timer()
	{
		m_timer_active = true;
		dispatch(this).source_set_event_handler <&status_logger::update_mt>(*m_message_timer);
		dispatch_resume(*m_message_timer);
	}
	
	
	void status_logger::clear_line_mt()
	{
		if (m_need_clear_line)
		{
			m_need_clear_line = false;
			//std::cerr << '\r';
			std::cerr << "\33[2K\r";
		}
	}
	
	
	void status_logger::handle_window_size_change_mt()
	{
		struct winsize ws;
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws);
		m_window_width = ws.ws_col;
		if (!m_window_width)
			m_window_width = 80;
		update_mt();
	}
}
