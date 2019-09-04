/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <vcf2multialign/preprocess/variant_graph_partitioner.hh>
#include <vcf2multialign/utility/check_ploidy.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/variant_format.hh>
#include "find_optimal_cut_positions.hh"


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
}


namespace vcf2multialign {
	
	void find_optimal_cut_positions(
		char const *reference_path,
		char const *variants_path,
		char const *output_path,
		char const *reference_seq_name,
		char const *chr_name,
		std::vector <std::string> const &field_names_for_filter_if_set,
		std::size_t const minimum_subgraph_distance,
		bool const should_overwrite_files
	)
	{
		vector_type reference;
		lb::file_ostream output_positions_stream;
		
		// Open the files.
		{
			auto const mode(lb::make_writing_open_mode({
				lb::writing_open_mode::CREATE,
				(should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
			}));
			lb::open_file_for_writing(output_path, output_positions_stream, mode);
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
		logging_variant_processor_delegate delegate;
		variant_graph_partitioner partitioner(
			delegate,
			reader,
			reference,
			chr_name,
			donor_count,
			chr_count,
			minimum_subgraph_distance
		);
		variant_graph_partitioner::cut_position_list cut_positions;
		partitioner.partition(field_names_for_filter_if_set, cut_positions);
		
		// Output.
		cereal::PortableBinaryOutputArchive archive(output_positions_stream);
		archive(cut_positions);
	}
}
