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
		
		virtual void finish() {}
		virtual void handle_variant(variant &var) {};
	};
	
	
	class variant_handler : public variant_buffer_delegate, public variant_handler_delegate
	{
	protected:
		dispatch_ptr <dispatch_queue_t>					m_main_queue{};
		dispatch_ptr <dispatch_queue_t>					m_parsing_queue{};
		
		variant_handler_delegate						*m_delegate{this};
		error_logger									*m_error_logger{};
		
		vector_type	const								*m_reference{};
		
		variant_buffer									m_variant_buffer;
		variant_set const								*m_skipped_variants{};
		std::set <size_t>								m_valid_alts;
		
		sv_handling										m_sv_handling_method{};
		std::size_t										m_i{0};
		bool											m_check_alts{true};
		
	public:
		variant_handler(
			dispatch_ptr <dispatch_queue_t> const &main_queue,
			dispatch_ptr <dispatch_queue_t> const &parsing_queue,
			vcf_reader &vcf_reader_,
			vector_type const &reference,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants,
			error_logger &error_logger
		):
			m_main_queue(main_queue),
			m_parsing_queue(parsing_queue),
			m_error_logger(&error_logger),
			m_reference(&reference),
			m_variant_buffer(vcf_reader_, main_queue, *this),
			m_skipped_variants(&skipped_variants),
			m_sv_handling_method(sv_handling_method)
		{
		}
		
		variant_handler() = default;
		
	public:
		variant_buffer &get_variant_buffer() { return m_variant_buffer; }
		void set_delegate(variant_handler_delegate &delegate) { m_delegate = &delegate; }
		bool is_valid_alt(uint8_t const alt_idx) const { return 0 < m_valid_alts.count(alt_idx); }
		std::set <size_t> const &valid_alts() const { return m_valid_alts; }
		
		void process_variants();
		void enumerate_genotype(
			variant &var,
			std::size_t const sample_no,
			std::function <void(uint8_t, std::size_t, bool)> const &cb
		);
		

	protected:
		virtual void handle_variant(variant &var) override;
		virtual void finish() override;
		
		bool check_alt_seq(std::string const &alt) const;
		void fill_valid_alts(variant const &var);
		void set_check_alts(bool const should_check) { m_check_alts = should_check; }
	};
}

#endif
