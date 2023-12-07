/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <filesystem>
#include <libbio/fasta_reader.hh>
#include <vcf2multialign/output.hh>
#include <vcf2multialign/variant_graph.hh>

namespace fs	= std::filesystem;
namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef v2m::founder_sequence_greedy_output::ploidy_matrix	ploidy_matrix;
	
	
	struct build_variant_graph_delegate final : public v2m::build_graph_delegate
	{
		bool should_include(std::string_view const sample_name, v2m::variant_graph::ploidy_type const chrom_copy_idx) const override { return true; }
		
		void report_overlapping_alternative(
			std::uint64_t const lineno,
			v2m::variant_graph::position_type const ref_pos,
			std::vector <std::string_view> const &var_id,
			std::string_view const sample_name,
			v2m::variant_graph::ploidy_type const chrom_copy_idx,
			std::uint32_t const gt
		) override
		{
			std::cerr << "Overlapping alternative alleles. Sample: " << sample_name << " chromosome copy: " << chrom_copy_idx << " current variant position: " << ref_pos << " genotype: " << gt << '\n';
		}
		
		bool ref_column_mismatch(std::uint64_t const var_idx, position_type const pos, std::string_view const expected, std::string_view const actual) override
		{
			FAIL("REF column contents do not match the reference sequence in variant " << var_idx << ", position " << pos << ". Expected: “" << expected << "” Actual: “" << actual << "”");
			return false;
		}
	};
	
	
	struct output_delegate final : public v2m::output_delegate
	{
		void will_handle_sample(std::string const &sample, sample_type const sample_idx, ploidy_type const chr_copy_idx) override {}
		void will_handle_founder_sequence(sample_type const idx) override {}
		void handled_sequences(sequence_count_type const sequence_count) override {}
		void exit_subprocess(v2m::subprocess_type &proc) override {}
		void handled_node(v2m::variant_graph::node_type const node) override {}
	};
	
	
	void test_founders(
		char const * const vcf_name,
		char const * const fasta_name,
		v2m::cut_position_vector const &expected_cut_positions,
		ploidy_matrix const &expected_matchings,
		char const * const expected_output
	)
	{
		INFO("VCF: " << vcf_name);
		INFO("FASTA: " << fasta_name);
		
		fs::path const base_path("test-files/founder-sequences");
		
		v2m::sequence_type ref_seq;
		REQUIRE(lb::read_single_fasta_sequence(base_path / fasta_name, ref_seq, nullptr));
		
		v2m::variant_graph graph;
		
		{
			build_variant_graph_delegate delegate;
			v2m::build_graph_statistics stats;
			v2m::build_variant_graph(ref_seq, base_path / vcf_name, "1", graph, stats, delegate);
		}
		
		{
			output_delegate delegate;
			v2m::founder_sequence_greedy_output output(nullptr, nullptr, true, false, false, delegate);
			REQUIRE(output.find_cut_positions(graph, 0));
			REQUIRE(expected_cut_positions == output.cut_positions());
			REQUIRE(output.find_matchings(graph, 2));
			REQUIRE(expected_matchings == output.assigned_samples());
			
			std::stringstream os;
			output.output_a2m(ref_seq, graph, os);
			REQUIRE(expected_output == os.view());
		}
	}
}


SCENARIO("Founder sequences can be generated from predef test input")
{
	GIVEN("A VCF file and a reference (1, 1)")
	{
		auto const * const expected_output(
			">REF\n"
			"CAA-AACTT-CCCGG-\n"
			">1\n"
			"AAA-AACTT-CCAGG-\n"
			">2\n"
			"CAA-AATTT-CCTGG-\n"
		);
		test_founders("test-1.vcf", "test-1.fa", {0, 1, 3, 5}, {{0, 6, 6, 3, 5, 8}, 3}, expected_output);
	}
	
	
	GIVEN("A VCF file and a reference (1, 1-2)")
	{
		auto const * const expected_output(
			">REF\n"
			"CAA-AACTT-CCCGG-AAAA\n"
			">1\n"
			"AAA-AACTT-CCAGG-AAAA\n"
			">2\n"
			"CAA-AATTT-CCTGG-AAAA\n"
		);
		test_founders("test-1.vcf", "test-1-2.fa", {0, 1, 3, 6}, {{0, 6, 6, 3, 5, 8}, 3}, expected_output);
	}
	
	
	GIVEN("A VCF file and a reference (2, 2)")
	{
		auto const * const expected_output(
			">REF\n"
			"CAA-CTTCG-G\n"
			">1\n"
			"CAA-CTTGG-G\n"
			">2\n"
			"AAA-CTGGGGG\n"
		);
		test_founders("test-2.vcf", "test-2.fa", {0, 3, 5}, {{6, 8, 0, 7}, 2}, expected_output);
	}
	
	
	GIVEN("A VCF file and a reference (3, 3)")
	{
		auto const * const expected_output(
			">REF\n"
			"CAA-CTT-CGG-\n"
			">1\n"
			"AAA-CTT-AGG-\n"
			">2\n"
			"CAA-TTT-TGG-\n"
		);
		test_founders("test-3.vcf", "test-3.fa", {0, 1, 2, 3}, {{0, 6, 6, 3, 5, 8}, 3}, expected_output);
	}
	
	
	GIVEN("A VCF file and a reference (4, 4)")
	{
		auto const * const expected_output(
			">REF\n"
			"TTTCAA-AACTT-CCCGG-\n"
			">1\n"
			"TTTAAA-AACTT-CCAGG-\n"
			">2\n"
			"TTTCAA-AATTT-CCTGG-\n"
		);
		test_founders("test-4.vcf", "test-4.fa", {0, 2, 4, 6}, {{0, 6, 6, 3, 5, 8}, 3}, expected_output);
	}
}
