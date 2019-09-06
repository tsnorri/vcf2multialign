/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCF_PROCESSOR_HH
#define VCF2MULTIALIGN_VCF_PROCESSOR_HH

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {

	class vcf_processor
	{
	protected:
		libbio::mmap_handle <char>	m_vcf_handle;
		libbio::vcf_mmap_input		m_vcf_input;
		libbio::vcf_reader			m_vcf_reader;
		
		
	public:
		vcf_processor():
			m_vcf_input(m_vcf_handle),
			m_vcf_reader(m_vcf_input)
		{
		}
		
		libbio::vcf_reader &vcf_reader() { return m_vcf_reader; }
		void open_variants_file(char const *variant_file_path);
		void prepare_reader();
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
