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
#include <vcf2multialign/preprocessing_result.hh>
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
	
	
	class variant_handler : public variant_buffer_delegate
	{
	protected:
		class variant_buffer_container
		{
		protected:
			variant_handler	*m_handler{nullptr};
			variant_buffer	m_variant_buffer;
			
		public:
			class variant_buffer const &variant_buffer() const { return m_variant_buffer; }
			class variant_buffer &variant_buffer() { return m_variant_buffer; }
			
			variant_buffer_container() = default;
			
			variant_buffer_container(variant_handler &handler): m_handler(&handler) {}
			
			// Use perfect forwarding for the remaining arguments.
			template <typename ... t_args>
			variant_buffer_container(variant_handler &handler, t_args && ... args):
				m_handler(&handler),
				m_variant_buffer(std::forward <t_args> (args)...)
			{
				m_variant_buffer.set_delegate(*m_handler);
			}
			
			variant_buffer_container(variant_buffer_container const &other):
				m_variant_buffer(other.m_variant_buffer)
			{
				m_variant_buffer.set_delegate(*m_handler);
			}
			
			variant_buffer_container(variant_buffer_container &&other):
				m_variant_buffer(std::move(other.m_variant_buffer))
			{
				m_variant_buffer.set_delegate(*m_handler);
			}
			
			variant_buffer_container &operator=(variant_buffer_container const &other) &
			{
				m_variant_buffer = other.m_variant_buffer;
				m_variant_buffer.set_delegate(*m_handler);
				return *this;
			}
			
			variant_buffer_container &operator=(variant_buffer_container &&other) &
			{
				m_variant_buffer = std::move(other.m_variant_buffer);
				m_variant_buffer.set_delegate(*m_handler);
				return *this;
			}
		};
		
	protected:
		dispatch_ptr <dispatch_queue_t>					m_parsing_queue{};
		
		variant_handler_delegate						*m_delegate{};
		error_logger									*m_error_logger{};
		preprocessing_result const						*m_preprocessing_result{};

		variant_buffer_container						m_vbc{*this};
		
		sv_handling										m_sv_handling_method{};
		
	public:
		variant_handler(
			error_logger &error_logger,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,	// Needs to be serial.
			dispatch_ptr <dispatch_queue_t> const &parsing_queue,	// May be concurrent since only variant_buffer's read_input is called there.
			class vcf_reader &vcf_reader,
			preprocessing_result const &result,
			sv_handling const sv_handling_method
		):
			m_parsing_queue(parsing_queue),
			m_error_logger(&error_logger),
			m_preprocessing_result(&result),
			m_vbc(*this, vcf_reader, worker_queue, *this),
			m_sv_handling_method(sv_handling_method)
		{
		}
		
		variant_handler() = default;
		virtual ~variant_handler() {}
		
		class variant_buffer const &variant_buffer() const { return m_vbc.variant_buffer(); }
		class variant_buffer &variant_buffer() { return m_vbc.variant_buffer(); }
		void set_delegate(variant_handler_delegate &delegate) { m_delegate = &delegate; }
		
		void process_variants();
		void enumerate_sample_genotypes(
			variant const &var,
			std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
		);

	protected:
		// variant_buffer_delegate.
		virtual void handle_variant(variant &var) override;
		virtual void finish() override;
	};
}

#endif
