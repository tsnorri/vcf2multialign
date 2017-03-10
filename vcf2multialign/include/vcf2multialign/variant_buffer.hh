/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_BUFFER_HH
#define VCF2MULTIALIGN_VARIANT_BUFFER_HH

#include <string>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {

	struct buffered_variant
	{
		variant		var;
		std::string	line;

		buffered_variant(std::size_t const sample_count = 0):
			var(sample_count)
		{
		}


		buffered_variant(std::size_t const sample_count, std::set <std::string> const &requested_format_fields):
			var(sample_count, requested_format_fields)
		{
		}
	};


	class variant_buffer
	{
	public:
		typedef std::vector <buffered_variant> variant_list;

	protected:
		vcf_reader			*m_reader{};
		variant_list		m_variant_list;
		buffered_variant	m_model_variant;
		std::size_t			m_list_ptr{0};

	public:
		variant_buffer():
			m_reader(nullptr),
			m_model_variant(0)
		{
		}
		
		variant_buffer(vcf_reader &reader, std::size_t const sample_count, std::set <std::string> const &requested_format_fields):
			m_reader(&reader),
			m_model_variant(sample_count, requested_format_fields)
		{
		}

		variant_buffer(vcf_reader &reader, std::size_t const sample_count = 0):
			m_reader(&reader),
			m_model_variant(sample_count)
		{
		}

		void add_format_field(std::string const &name) { m_model_variant.var.add_format_field(name); }
		std::pair <variant_list::iterator, variant_list::iterator> variant_range();
		void fill_buffer();
	};
}

#endif
