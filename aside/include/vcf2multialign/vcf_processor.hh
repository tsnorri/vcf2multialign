/*
 * Copyright (c) 2019-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCF_PROCESSOR_HH
#define VCF2MULTIALIGN_VCF_PROCESSOR_HH

#include <libbio/vcf/vcf_reader.hh>
#include <memory>
#include <vcf2multialign/types.hh>


namespace vcf2multialign::detail {
	
	struct vcf_input
	{
		virtual ~vcf_input() {}
	};
}


namespace vcf2multialign {

	class vcf_processor
	{
	protected:
		typedef std::unique_ptr <detail::vcf_input>					input_ptr;
		typedef std::unique_ptr <libbio::vcf::variant_validator>	variant_validator_ptr;
		
	protected:
		input_ptr				m_vcf_input;
		variant_validator_ptr	m_variant_validator;
		libbio::vcf::reader		m_vcf_reader;
		
	public:
		libbio::vcf::reader &vcf_reader() { return m_vcf_reader; }
		void open_variants_file(char const *variant_file_path, char const *bed_file_path = nullptr);
		void prepare_reader();
		void setup_empty_input();
	};
	
	
	class reference_controller
	{
	protected:
		vector_type					m_reference;
		
	public:
		void read_reference(char const *reference_path, char const *reference_seq_name);
	};
	
	
	class output_stream_controller
	{
	protected:
		libbio::file_ostream		m_output_stream;
		
	public:
		output_stream_controller() = default;
		void open_output_file(char const *output_path, bool const should_overwrite_files);
	};
}

#endif
