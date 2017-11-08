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
	
	struct status_logger_delegate
	{
		virtual ~status_logger_delegate() {}
		virtual std::size_t step_count() const = 0;
		virtual std::size_t current_step() const = 0;
	};
	
	
	class status_logger
	{
	public:
		typedef std::chrono::time_point <std::chrono::steady_clock> time_point_type;
		
	protected:
		enum indicator_type : uint8_t
		{
			none = 0,
			counter,
			progress_bar
		};
		
	protected:
		status_logger_delegate				*m_delegate{nullptr};
		dispatch_ptr <dispatch_source_t>	m_message_timer{nullptr};
		dispatch_ptr <dispatch_source_t>	m_signal_source{nullptr};
		std::mutex							m_message_mutex;
		std::string							m_message;
		time_point_type						m_start_time{};
		std::size_t							m_window_width{0};
		std::size_t							m_message_length{0};
		std::atomic <indicator_type>		m_indicator_type{none};
		std::atomic_bool					m_timer_active{false};
		bool								m_need_clear_line{false};
		
	public:
		status_logger() = default;
		
		status_logger(status_logger_delegate &delegate):
			m_delegate(&delegate)
		{
		}
		
		void set_delegate(status_logger_delegate &delegate) { m_delegate = &delegate; }
		void set_message(std::string const &message);
		void install();
		void uninstall();
		void finish_logging();
		void log_message_counting(std::string const &message);
		void log_message_progress_bar(std::string const &message);
		void update();
		
		template <typename t_fn>
		void log(t_fn fn, bool clear_line = true)
		{
			dispatch_async_fn(dispatch_get_main_queue(), [this, fn{std::move(fn)}, clear_line](){
				if (clear_line)
					clear_line_mt();
				fn();
				update_mt();
			});
		}
		
	protected:
		void resume_timer();
		void update_mt();
		void clear_line_mt();
		void finish_logging_mt();
		void handle_window_size_change_mt();
	};
}

#endif
