/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/variant_format.hh>
#include <vcf2multialign/vcf_processor.hh>

namespace lb	= libbio;
namespace vcf	= libbio::vcf;


namespace {
	
	template <typename t_map, typename t_key, typename t_ptr_value>
	void add_to_map(t_map &map, t_key const &key, t_ptr_value const ptr_value)
	{
		map.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(key),
			std::forward_as_tuple(ptr_value)
		);
	}
}


namespace vcf2multialign {

	void reference_controller::read_reference(char const *reference_path, char const *reference_seq_name)
	{
		// Read the input FASTA.
		lb::mmap_handle <char> ref_handle;
		ref_handle.open(reference_path);
		read_single_fasta_seq(ref_handle, m_reference, reference_seq_name);
	}
	
	
	void vcf_processor::open_variants_file(char const *variant_file_path)
	{
		m_vcf_handle.open(variant_file_path);
		m_vcf_input.reset_range();
	}
	
	
	void output_stream_controller::open_output_file(char const *output_path, bool const should_overwrite_files)
	{
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		lb::open_file_for_writing(output_path, m_output_stream, mode);
	}
	
	
	void vcf_processor::prepare_reader()
	{
		vcf::add_reserved_genotype_keys(m_vcf_reader.genotype_fields());
	
		{
			auto &info_fields(m_vcf_reader.info_fields());
			vcf::add_reserved_info_keys(info_fields);
	
			// Add some fields used in 1000G.
			add_to_map(info_fields, "CIPOS",	new vcf_info_field_cipos());
			add_to_map(info_fields, "CIEND",	new vcf_info_field_ciend());
			add_to_map(info_fields, "SVLEN",	new vcf_info_field_svlen());
			add_to_map(info_fields, "SVTYPE",	new vcf_info_field_svtype());
		}
	
		m_vcf_reader.set_variant_format(new variant_format());
		m_vcf_reader.fill_buffer();
		m_vcf_reader.read_header();
	}
}
