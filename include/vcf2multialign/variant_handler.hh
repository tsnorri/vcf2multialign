/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_HANDLER_HH
#define VCF2MULTIALIGN_VARIANT_HANDLER_HH

#include <dispatch/dispatch.h>
#include <map>
#include <stack>
#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_buffer.hh>
#include <vcf2multialign/vcf_reader.hh>
#include <vcf2multialign/vector_source.hh>


namespace vcf2multialign {
	
	struct variant_handler_delegate
	{
		virtual ~variant_handler_delegate() {}
		
		virtual void prepare(vcf_reader &reader)
		{
			reader.set_parsed_fields(vcf_field::ALL);
		}
		
		virtual void handle_variant(variant &var) {};
		virtual void finish() {}
		
		virtual std::vector <uint8_t> const &valid_alts(std::size_t const lineno) const = 0;
		virtual bool is_valid_alt(std::size_t const lineno, uint8_t const alt_idx) const = 0;
	};
	
	
	class variant_handler_base : public variant_buffer_delegate
	{
	protected:
		dispatch_ptr <dispatch_queue_t>					m_parsing_queue{};
		
		variant_handler_delegate						*m_delegate{};
		error_logger									*m_error_logger{};
		
		vector_type	const								*m_reference{};
		
		variant_buffer									m_variant_buffer;
		variant_set const								*m_skipped_variants{};
		
		sv_handling										m_sv_handling_method{};
		
	public:
		variant_handler_base(
			dispatch_ptr <dispatch_queue_t> const &worker_queue,	// Needs to be serial.
			dispatch_ptr <dispatch_queue_t> const &parsing_queue,	// May be concurrent since only variant_buffer's read_input is called there.
			class vcf_reader &vcf_reader,
			vector_type const &reference,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants,
			error_logger &error_logger
		):
			m_parsing_queue(parsing_queue),
			m_error_logger(&error_logger),
			m_reference(&reference),
			m_variant_buffer(vcf_reader, worker_queue, *this),
			m_skipped_variants(&skipped_variants),
			m_sv_handling_method(sv_handling_method)
		{
		}
		
		variant_handler_base() = default;
		virtual ~variant_handler_base() {}
		
	protected:
		virtual void finish() = 0;
		virtual void handle_variant(variant &var) = 0;
	};
	
	
	class variant_handler : public variant_handler_base
	{
	protected:
		void finish_copy()
		{
			m_variant_buffer.set_delegate(*this);
		}
		
	public:
		using variant_handler_base::variant_handler_base;
		
		variant_handler(variant_handler const &other):
			variant_handler_base(other)
		{
			finish_copy();
		}
		
		variant_handler &operator=(variant_handler const &other) &
		{
			variant_handler_base::operator=(other);
			finish_copy();
			return *this;
		}
		
	public:
		class variant_buffer &variant_buffer() { return m_variant_buffer; }
		void set_delegate(variant_handler_delegate &delegate) { m_delegate = &delegate; }
		
		void process_variants();
		void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		);

	protected:
		virtual void handle_variant(variant &var) override;
		virtual void finish() override;
	};
}

#endif
