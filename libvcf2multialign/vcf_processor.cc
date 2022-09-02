/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <vcf2multialign/bed_reader_delegate.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/variant_format.hh>
#include <vcf2multialign/vcf_processor.hh>

namespace io	= boost::iostreams;
namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace {
	
	template <std::size_t t_size>
	inline bool ends_with(std::string_view const &sv, char const (&ending_arr)[t_size])
	{
		auto const sv_size(sv.size());
		if (sv_size < t_size - 1) // Remove the trailing nul character.
			return false;

		std::string_view const ending(ending_arr, t_size - 1);
		auto const suffix(sv.substr(sv_size - t_size + 1));
		return (suffix == ending);
	}
	
	
	struct vcf_empty_input final : public v2m::detail::vcf_input
	{
		vcf::empty_input input{};
	};
	
	
	struct vcf_mmap_input final : public v2m::detail::vcf_input
	{
		vcf::mmap_input	input{};
	};
	
	
	// FIXME: make libbio handle compressed input or read input from a pipe.
	struct vcf_compressed_stream_input final : public v2m::detail::vcf_input
	{
		typedef vcf::stream_input <
			io::filtering_stream <io::input>
		> filtering_stream_input;
		
		filtering_stream_input	input;
		lb::file_istream		compressed_input_stream;
	};
	
	
	void output_subfield_description(std::ostream &os, vcf::subfield_base const &field)
	{
		vcf::output_vcf_value(os, field.metadata_value_type());
		vcf::output_vcf_value(os, field.number());
		if (field.value_type_is_vector())
			os << 'V';
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
	
	
	void vcf_processor::open_variants_file(char const *variant_file_path, char const *bed_file_path)
	{
		std::string_view const sv(variant_file_path);
		if (ends_with(sv, ".gz")) // std::string_view::ends_with is a C++20 addition.
		{
			auto ptr(std::make_unique <vcf_compressed_stream_input>());
			lb::open_file_for_reading(variant_file_path, ptr->compressed_input_stream);
			
			auto &filtering_stream(ptr->input.stream());
			filtering_stream.push(boost::iostreams::gzip_decompressor());
			filtering_stream.push(ptr->compressed_input_stream);
			filtering_stream.exceptions(std::istream::badbit);
			
			m_vcf_reader.set_input(ptr->input);
			m_vcf_input = std::move(ptr);
		}
		else
		{
			auto ptr(std::make_unique <vcf_mmap_input>());
			ptr->input.handle().open(variant_file_path);
			
			m_vcf_reader.set_input(ptr->input);
			m_vcf_input = std::move(ptr);
		}
		
		if (bed_file_path)
		{
			auto validator(std::make_unique <vcf::region_variant_validator>(true));
			lb::bed_reader bed_reader;
			bed_reader_delegate delegate(validator->regions());
			bed_reader.read_regions(bed_file_path, delegate);
			m_variant_validator = std::move(validator);
		}
		else
		{
			m_variant_validator.reset(new vcf::region_variant_validator(false));
		}
	}
	
	
	void vcf_processor::setup_empty_input()
	{
		auto ptr(std::make_unique <vcf_empty_input>());
		m_vcf_reader.set_input(ptr->input);
		m_vcf_input = std::move(ptr);
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
			vcf::add_subfield <vcf_info_field_cipos>  (info_fields, "CIPOS");
			vcf::add_subfield <vcf_info_field_ciend>  (info_fields, "CIEND");
			vcf::add_subfield <vcf_info_field_svlen>  (info_fields, "SVLEN");
			vcf::add_subfield <vcf_info_field_svtype> (info_fields, "SVTYPE");
		}
		
		m_vcf_reader.set_variant_format(new variant_format());
		
		m_vcf_reader.read_header();
	}
}
