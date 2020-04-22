/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_PROCESSOR_HH
#define VCF2MULTIALIGN_VARIANT_PROCESSOR_HH

#include <libbio/vcf/subfield.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_processor_delegate.hh>


namespace vcf2multialign {

	class variant_processor
	{
	public:
		enum class variant_check_status : std::uint8_t {
			PASS,
			ERROR,
			FATAL_ERROR
		};
		
		typedef std::vector <libbio::vcf::info_field_base *> vcf_info_field_vector;
		
	protected:
		libbio::vcf::reader								*m_reader{};
		vector_type const								*m_reference{};
		std::string const								*m_chromosome_name{};
		
	public:
		variant_processor() = default;
		
		variant_processor(
			libbio::vcf::reader &reader,
			vector_type const &reference,
			std::string const &chromosome_name
		):
			m_reader(&reader),
			m_reference(&reference),
			m_chromosome_name(&chromosome_name)
		{
		}
		
		virtual ~variant_processor() {}
		
		libbio::vcf::reader &reader() { return *m_reader; }
		vector_type const &reference() { return *m_reference; }
		
		void fill_filter_by_assigned(
			std::vector <std::string> const &field_names_for_filter_by_assigned,
			vcf_info_field_vector &filter_by_assigned,
			variant_processor_delegate &delegate
		);
		
		variant_check_status check_variant(
			libbio::vcf::transient_variant const &var,
			vcf_info_field_vector const &filter_by_assigned,
			variant_processor_delegate &delegate
		);
	};
}

#endif
