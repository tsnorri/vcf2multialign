/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <vcf2multialign/preprocess/variant_preprocessor.hh>
#include <vcf2multialign/utility/check_ploidy.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/variant_format.hh>
#include "preprocess_variants.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	std::tuple <std::size_t, std::uint8_t> check_ploidy(lb::vcf_reader &vcf_reader)
	{
		v2m::ploidy_map ploidy_map;
		v2m::check_ploidy(vcf_reader, ploidy_map);
		
		auto const donor_count(ploidy_map.size());
		auto const ploidy(ploidy_map.begin()->second);
		
		bool can_continue(true);
		for (auto const [sample_no, sample_ploidy] : ploidy_map)
		{
			if (sample_ploidy != ploidy)
			{
				std::cerr << "ERROR: Ploidy for sample " << sample_no << " was " << sample_ploidy << ", expected " << ploidy << ".\n";
				can_continue = false;
			}
		}
		
		if (!can_continue)
			libbio_fail("Varying ploidy is not currently supported.");
		
		return std::tuple <std::size_t, std::uint8_t>(donor_count, ploidy);
	}
	
	
	template <typename t_map, typename t_key, typename t_ptr_value>
	void add_to_map(t_map &map, t_key const &key, t_ptr_value const ptr_value)
	{
		map.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(key),
			std::forward_as_tuple(ptr_value)
		);
	}
	
	
	struct graph_output_delegate final : public v2m::variant_preprocessor_delegate
	{
		std::size_t subgraph_number{};
		
		virtual void variant_preprocessor_no_field_for_identifier(std::string const &identifier) override
		{
			std::cerr << "WARNING: Did not find a field for identifier “" << identifier << "”.\n";
		}
		
		
		virtual void variant_preprocessor_found_variant_with_position_greater_than_reference_length(lb::transient_variant const &var) override
		{
			std::cerr << "ERROR: Found a variant with a position greater than the reference length on line " << var.lineno() << "\n";
		}
		
		
		virtual void variant_preprocessor_found_variant_with_no_suitable_alts(lb::transient_variant const &var) override
		{
			std::cerr << "Line " << var.lineno() << ": Variant has no ALTs that could be handled.\n";
		}
		
		
		virtual void variant_preprocessor_found_filtered_variant(lb::transient_variant const &var, lb::vcf_info_field_base const &field) override
		{
			std::cerr << "Line " << var.lineno() << ": Variant has the field '" << field.get_metadata()->get_id() << "' set; skipping.\n";
		}
		
		
		virtual void variant_preprocessor_found_variant_with_ref_mismatch(lb::transient_variant const &var, std::string_view const &ref_sub) override
		{
			std::cerr << "WARNING: reference column mismatch on line " << var.lineno() << ": expected '" << ref_sub << "', got '" << var.ref() << "'\n";
		}
		
		
		virtual void variant_preprocessor_will_handle_subgraph(lb::variant const &first_var, std::size_t const variant_count, std::size_t const path_count) override
		{
			std::cout << subgraph_number++ << '\t' << first_var.lineno() << '\t' << variant_count << '\t' << path_count << '\n';
		}
	};
}


namespace vcf2multialign {
	
	void preprocess_variants(
		char const *reference_path,
		char const *variants_path,
		char const *output_variants_path,
		char const *reference_seq_name,
		char const *chr_name,
		std::vector <std::string> const &field_names_for_filter_if_set,
		std::size_t const minimum_subgraph_distance,
		bool const should_overwrite_files
	)
	{
		vector_type reference;
		lb::file_ostream output_graph_stream;
		
		// Open the files.
		{
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			lb::open_file_for_writing(output_variants_path, output_graph_stream, mode);
		}
		
		lb::mmap_handle <char> vcf_handle;
		vcf_handle.open(variants_path);
		
		lb::mmap_handle <char> ref_handle;
		ref_handle.open(reference_path);
		
		// Read the input FASTA.
		read_single_fasta_seq(ref_handle, reference, reference_seq_name);

		// Create a VCF reader.
		lb::vcf_mmap_input vcf_input(vcf_handle);
		lb::vcf_reader reader(vcf_input);
		lb::add_reserved_genotype_keys(reader.genotype_fields());
		
		{
			auto &info_fields(reader.info_fields());
			lb::add_reserved_info_keys(info_fields);
			
			// Add some fields used in 1000G.
			add_to_map(info_fields, "CIPOS",	new vcf_info_field_cipos());
			add_to_map(info_fields, "CIEND",	new vcf_info_field_ciend());
			add_to_map(info_fields, "SVLEN",	new vcf_info_field_svlen());
			add_to_map(info_fields, "SVTYPE",	new vcf_info_field_svtype());
		}
		
		reader.set_variant_format(new variant_format());
		
		// Parse.
		reader.set_parsed_fields(lb::vcf_field::ALL);
		reader.set_input(vcf_input);
		reader.fill_buffer();
		reader.read_header();
		
		// Check the ploidy.
		ploidy_map ploidy;
		check_ploidy(reader, ploidy);
		auto const donor_count(ploidy.size());
		if (!donor_count)
		{
			std::cerr << "WARNING: No donors found." << std::endl;
			return;
		}
		auto const chr_count(ploidy.begin()->second);
		
		// Process the variants.
		graph_output_delegate delegate;
		std::cout << "SUBGRAPH_NUMBER" << '\t' << "LINE_NUMBER" << '\t' << "VARIANT_COUNT" << '\t' << "PATH_COUNT" << '\n';
		variant_preprocessor processor(delegate, reader, reference, chr_name, donor_count, chr_count, minimum_subgraph_distance);
		processor.process(field_names_for_filter_if_set);
		
		// Output.
		cereal::PortableBinaryOutputArchive archive(output_graph_stream);
		archive(processor.variant_graph());
	}
}
