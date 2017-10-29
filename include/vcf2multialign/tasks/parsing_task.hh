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
	class parsing_task_base : public task
	{
	protected:
		status_logger						*m_status_logger{nullptr};
		error_logger						*m_error_logger{nullptr};
		vcf_reader							m_vcf_reader;
		
	public:
		parsing_task_base() = default;
		
		parsing_task_base(
			status_logger &status_logger,
			error_logger &error_logger,
			class vcf_reader const &vcf_reader
		):
			m_status_logger(&status_logger),
			m_error_logger(&error_logger),
			m_vcf_reader(vcf_reader)
		{
		}
		
		parsing_task_base(
			status_logger &status_logger,
			error_logger &error_logger,
			class vcf_reader &&vcf_reader
		):
			m_status_logger(&status_logger),
			m_error_logger(&error_logger),
			m_vcf_reader(std::move(vcf_reader))
		{
		}
		
		virtual ~parsing_task_base() {}
		virtual void execute() = 0;
		
		class vcf_reader &vcf_reader() { return m_vcf_reader; }
	};
	
	
	class parsing_task : public parsing_task_base
	{
	protected:
		void finish_copy()
		{
			m_vcf_reader.reset();
		}
		
	public:
		using parsing_task_base::parsing_task_base;
		
		parsing_task(parsing_task const &other):
			parsing_task_base(other)
		{
			finish_copy();
		}
		
		parsing_task(parsing_task &&other):
			parsing_task_base(std::move(other))
		{
		}
		
		parsing_task &operator=(parsing_task const &other) &
		{
			parsing_task_base::operator=(other);
			finish_copy();
			return *this;
		}
		
		parsing_task &operator=(parsing_task &&other) &
		{
			parsing_task_base::operator=(std::move(other));
			return *this;
		}
	};
	
	
	class parsing_task_vh_base : public parsing_task, public variant_handler_delegate
	{
	protected:
		alt_checker const	*m_checker{nullptr};
		variant_handler		m_variant_handler;
	
	public:
		parsing_task_vh_base() = default;
		
		parsing_task_vh_base(
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
			m_checker(&checker),
			m_variant_handler(
				worker_queue,
				dispatch_ptr <dispatch_queue_t>(dispatch_get_main_queue()),
				m_vcf_reader,
				reference,
				sv_handling_method,
				skipped_variants,
				error_logger
			)
		{
			m_variant_handler.set_delegate(*this);
		}
	};
	
	
	class parsing_task_vh : public parsing_task_vh_base
	{
	protected:
		void finish_copy_or_move()
		{
			m_variant_handler.variant_buffer().set_vcf_reader(m_vcf_reader);
			m_variant_handler.set_delegate(*this);
		}
	
	public:
		using parsing_task_vh_base::parsing_task_vh_base;
		
		parsing_task_vh(parsing_task_vh const &other):
			parsing_task_vh_base(other)
		{
			finish_copy_or_move();
		}
		
		parsing_task_vh(parsing_task_vh &&other):
			parsing_task_vh_base(std::move(other))
		{
			finish_copy_or_move();
		}
		
		parsing_task_vh &operator=(parsing_task_vh const &other) &
		{
			parsing_task_vh_base::operator=(other);
			finish_copy_or_move();
			return *this;
		}
		
		parsing_task_vh &operator=(parsing_task_vh &&other) &
		{
			parsing_task_vh_base::operator=(std::move(other));
			finish_copy_or_move();
			return *this;
		}
		
		// variant_handler_delegate
		virtual std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const override { return m_checker->valid_alts(lineno); }
		virtual bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const override { return m_checker->is_valid_alt(lineno, alt_idx); }
	};
}

#endif
