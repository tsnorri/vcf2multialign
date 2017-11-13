/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_PARSING_TASK_HH
#define VCF2MULTIALIGN_TASKS_PARSING_TASK_HH

#include <vcf2multialign/alt_checker.hh>
#include <vcf2multialign/file_handling.hh>
#include <vcf2multialign/status_logger.hh>
#include <vcf2multialign/tasks/task.hh>
#include <vcf2multialign/variant_handler.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {
	
	// Since vcf_reader stores a pointer to the input stream,
	// copying and moving will leave a dangling pointer.
	// Use a base class to create default implementations of
	// copy and move constructors and copy and move assignment
	// operators and set the pointer in a subclass.
	//
	// variant_handler has the following pointers.
	// At the monent, the ones marked need to be handled by
	// the owner of parsing_task.
	//  [x] variant_handler_delegate	*m_delegate;
	//  [x] error_logger				*m_error_logger;
	//      vector_type	const			*m_reference;
	//  [x] variant_set const			*m_skipped_variants;
	//
	// variant_handler has a variant_buffer, which has the following pointers (in a struct):
	//      vcf_reader					*m_reader;
	//      variant_buffer_delegate		*m_delegate;
	//
	// vcf_reader has the following pointer:
	//      vcf_input				*m_input; (owned by generate_context)
	class parsing_task : public task
	{
	protected:
		struct vcf_reader_container
		{
			class vcf_reader vcf_reader;
			
			vcf_reader_container() = default;
			vcf_reader_container(class vcf_reader const &vcf_reader_): vcf_reader(vcf_reader_) { vcf_reader.reset(); }
			
			vcf_reader_container(vcf_reader_container const &other):
				vcf_reader(other.vcf_reader)
			{
				vcf_reader.reset();
			}
			
			vcf_reader_container(vcf_reader_container &&other):
				vcf_reader(std::move(other.vcf_reader))
			{
				vcf_reader.reset();
			}
			
			vcf_reader_container &operator=(vcf_reader_container const &other) &
			{
				vcf_reader = other.vcf_reader;
				vcf_reader.reset();
				return *this;
			}
			
			vcf_reader_container &operator=(vcf_reader_container &&other) &
			{
				vcf_reader = std::move(other.vcf_reader);
				vcf_reader.reset();
				return *this;
			}
		};
		
	protected:
		status_logger						*m_status_logger{nullptr};
		error_logger						*m_error_logger{nullptr};
		vcf_reader_container				m_vrc;
		
	public:
		parsing_task() = default;
		
		parsing_task(
			status_logger &status_logger,
			error_logger &error_logger,
			class vcf_reader const &vcf_reader
		):
			m_status_logger(&status_logger),
			m_error_logger(&error_logger),
			m_vrc(vcf_reader)
		{
		}
		
		parsing_task(
			status_logger &status_logger,
			error_logger &error_logger,
			class vcf_reader &&vcf_reader
		):
			m_status_logger(&status_logger),
			m_error_logger(&error_logger),
			m_vrc(std::move(vcf_reader))
		{
		}
		
		virtual ~parsing_task() {}
		virtual void execute() = 0;
		
		class vcf_reader const &vcf_reader() const { return m_vrc.vcf_reader; }
		class vcf_reader &vcf_reader() { return m_vrc.vcf_reader; }
	};
	
	
	class parsing_task_vh : public parsing_task, public variant_handler_delegate
	{
	protected:
		class variant_handler_container
		{
		protected:
			parsing_task_vh	*m_task{nullptr};
			variant_handler	m_variant_handler;
			
		public:
			class variant_handler const &variant_handler() const { return m_variant_handler; }
			class variant_handler &variant_handler() { return m_variant_handler; }
			
			variant_handler_container() = default;
			
			// Use perfect forwarding for the remaining arguments.
			template <typename ... t_args>
			variant_handler_container(parsing_task_vh &task, t_args && ... args):
				m_task(&task),
				m_variant_handler(std::forward <t_args> (args)...)
			{
			}
			
			variant_handler_container(variant_handler_container const &other):
				m_variant_handler(other.m_variant_handler)
			{
				finish_copy_or_move();
			}
			
			variant_handler_container(variant_handler_container &&other):
				m_variant_handler(std::move(other.m_variant_handler))
			{
				finish_copy_or_move();
			}
			
			variant_handler_container &operator=(variant_handler_container const &other) &
			{
				m_variant_handler = other.m_variant_handler;
				finish_copy_or_move();
				return *this;
			}
			
			variant_handler_container &operator=(variant_handler_container &&other) &
			{
				m_variant_handler = std::move(other.m_variant_handler);
				finish_copy_or_move();
				return *this;
			}
			
		protected:
			void finish_copy_or_move()
			{
				m_variant_handler.variant_buffer().set_vcf_reader(m_task->vcf_reader());
				m_variant_handler.set_delegate(*m_task);
			}
		};
		
	protected:
		alt_checker const			*m_alt_checker{nullptr};
		variant_handler_container	m_vhc;
	
	public:
		parsing_task_vh() = default;
		
		parsing_task_vh(
			dispatch_ptr <dispatch_queue_t> const &worker_queue,	// Needs to be serial.
			status_logger &status_logger,
			error_logger &error_logger,
			class vcf_reader const &vcf_reader,
			alt_checker const &checker,
			vector_type const &reference,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants
		):
			parsing_task(status_logger, error_logger, vcf_reader),
			m_alt_checker(&checker),
			m_vhc(
				*this,
				worker_queue,
				dispatch_ptr <dispatch_queue_t>(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0)),
				m_vrc.vcf_reader,
				reference,
				sv_handling_method,
				skipped_variants,
				error_logger
			)
		{
			m_vhc.variant_handler().set_delegate(*this);
		}
		
		virtual ~parsing_task_vh() {}
		class variant_handler const &variant_handler() const { return m_vhc.variant_handler(); }
		class variant_handler &variant_handler() { return m_vhc.variant_handler(); }
		
		// variant_handler_delegate
		virtual std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const override { return m_alt_checker->valid_alts(lineno); }
		virtual bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const override { return m_alt_checker->is_valid_alt(lineno, alt_idx); }
	};
}

#endif
