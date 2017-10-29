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
		virtual std::size_t record_count() const = 0;
		virtual std::size_t current_record() const = 0;
	};
	
	
	class status_logger
	{
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
		std::size_t							m_window_width{0};
		std::size_t							m_message_length{0};
		std::atomic <indicator_type>		m_indicator_type{none};
		bool								m_need_clear_line{false};
		
	public:
		status_logger() = default;
		
		status_logger(status_logger_delegate &delegate):
			m_delegate(&delegate)
		{
		}
		
		void set_delegate(status_logger_delegate &delegate) { m_delegate = &delegate; }
		void install();
		void uninstall();
		void finish_logging();
		void log_message_counting(std::string const &message);
		void log_message_progress_bar(std::string const &message);
		void update();
		
		template <typename t_fn>
		void log(t_fn fn)
		{
			dispatch_async_fn(dispatch_get_main_queue(), [this, fn{std::move(fn)}](){
				clear_line_mt();
				fn();
				update_mt();
			});
		}
		
	protected:
		void update_mt();
		void clear_line_mt();
		void finish_logging_mt();
		void handle_window_size_change_mt();
	};
}

#endif