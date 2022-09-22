/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <vcf2multialign/variant_graph/variant_graph.hh>
#include <vcf2multialign/variant_graph/variant_graph_generator.hh>
#include <vcf2multialign/utility/find_first_matching_variant.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/variant_format.hh>


namespace gen	= Catch::Generators;
namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;
namespace vgs	= vcf2multialign::variant_graphs;


namespace {
	
	struct variant_graph_generator_delegate final : public vgs::variant_graph_single_pass_generator_delegate
	{
		void variant_graph_generator_will_handle_subgraph(
			lb::vcf::variant const &first_var,
			std::size_t const variant_count,
			std::size_t const path_count
		) override
		{
		}
		
		void variant_processor_found_variant_with_chrom_id_mismatch(lb::vcf::transient_variant const &var) override {}
		void variant_processor_no_field_for_identifier(std::string const &identifier) override {}
		void variant_processor_found_variant_with_position_greater_than_reference_length(lb::vcf::transient_variant const &var) override {}
		void variant_processor_found_variant_with_no_suitable_alts(lb::vcf::transient_variant const &var) override {}
		void variant_processor_found_filtered_variant(lb::vcf::transient_variant const &var, lb::vcf::info_field_base const &field) override {}
		void variant_processor_found_variant_with_ref_mismatch(lb::vcf::transient_variant const &var, std::string_view const &ref_sub) override {}
		void variant_processor_found_matching_variant(libbio::vcf::transient_variant const &var) override {}
	};
	
	
	struct alt_edge
	{
		std::size_t	target_node{};
		std::string	label;
		
		alt_edge() = default;
		
		alt_edge(std::size_t const target_node_, std::string const &label_):
			target_node(target_node_),
			label(label_)
		{
		}
	};


	bool operator==(vgs::variant_graph::alt_edge const &lhs_edge, alt_edge const &rhs_edge)
	{
		return lhs_edge.target_node == rhs_edge.target_node && lhs_edge.label == rhs_edge.label;
	}
	
	
	struct node_description
	{
		std::vector <alt_edge>	alt_edges;
		std::string				ref;
		std::size_t				node{};
		std::size_t				pos{};
		std::size_t				aln_pos{};
		bool					begins_subgraph{};
		
		node_description() = default;
		
		node_description(
			std::size_t node_,
			std::size_t pos_,
			std::size_t aln_pos_,
			bool begins_subgraph_,
			std::string ref_,
			std::vector <alt_edge> alt_edges_
		):
			alt_edges(std::move(alt_edges_)),
			ref(std::move(ref_)),
			node(node_),
			pos(pos_),
			aln_pos(aln_pos_),
			begins_subgraph(begins_subgraph_)
		{
		}
	};
	
	
	class node_comparator
	{
	protected:
		std::vector <node_description>	m_node_descriptions;
		
	public:
		node_comparator() = default;
		
		node_comparator(std::vector <node_description> node_descriptions):
			m_node_descriptions(std::move(node_descriptions))
		{
		}
		
		void check_graph(vgs::variant_graph const &graph, v2m::vector_type const &reference)
		{
			// FIXME: path checks.
			
			vgs::variant_graph_walker walker(graph, reference);
			walker.setup();
			
			auto it(m_node_descriptions.begin());
			auto const end(m_node_descriptions.end());
			while (true)
			{
				bool current_node_begins_subgraph(false);
				auto const state(walker.advance_and_track_subgraph());
				switch (state)
				{
					case vgs::variant_graph_walker::state::NODE:
						break;
					
					case vgs::variant_graph_walker::state::SUBGRAPH_START_NODE:
						current_node_begins_subgraph = true;
						break;
					
					case vgs::variant_graph_walker::state::END:
						goto finish;
				}
				
				REQUIRE(it != end);
				
				auto const &node_desc(*it);
				CHECK(node_desc.node == walker.node());
				INFO("Node " << node_desc.node)
				CHECK(node_desc.pos == walker.ref_position());
				INFO("Ref pos " << node_desc.pos);
				CHECK(node_desc.aln_pos == walker.aligned_position());
				INFO("Aln pos " << node_desc.aln_pos);
				CHECK(node_desc.begins_subgraph == current_node_begins_subgraph);
				INFO("Begins subgraph " << node_desc.begins_subgraph);

				auto const &expected_alt_edges(node_desc.alt_edges);
				auto const &actual_alt_edges(walker.alt_edges());
				CHECK(expected_alt_edges.size() == actual_alt_edges.size());
				for (auto const &[edge_idx, tup] : rsv::enumerate(rsv::zip(expected_alt_edges, actual_alt_edges)))
				{
					auto const &[expected_edge, actual_edge] = tup;
					INFO("ALT edge " << edge_idx);
					INFO("Expected target: " << expected_edge.target_node);
					INFO("Actual target:   " << actual_edge.target_node);
					INFO("Expected label:  " << expected_edge.label);
					INFO("Actual label:    " << actual_edge.label);
					CHECK(expected_edge == actual_edge);
				}
				
				++it;
			}
			
		finish:
			CHECK(it == end);
		}
	};
	
	
	void test_variant_graph(char const *vcf_name, char const *fasta_name, node_comparator &cmp)
	{
		std::string const data_file_dir("test-files/variant-graph/");
		auto const vcf_path(data_file_dir + vcf_name);
		auto const fasta_path(data_file_dir + fasta_name);
		
		v2m::vector_type reference_seq;
		v2m::read_single_fasta_seq(fasta_path.c_str(), reference_seq, nullptr, false);
		
		lb::vcf::mmap_input vcf_input;
		vcf_input.handle().open(vcf_path);
		lb::vcf::reader vcf_reader(vcf_input);
		
		lb::vcf::add_reserved_info_keys(vcf_reader.info_fields());
		lb::vcf::add_reserved_genotype_keys(vcf_reader.genotype_fields());
		
		vcf_reader.set_variant_format(new v2m::variant_format());
		vcf_reader.read_header();
		
		// The parsed fields need to be set here so that the first matching variant gets parsed correctly.
		vcf_reader.set_parsed_fields(vcf::field::ALL);
		REQUIRE(v2m::find_first_matching_variant(vcf_reader, "1"));
		
		variant_graph_generator_delegate delegate;
		vgs::variant_graph_single_pass_generator generator(
			delegate,
			vcf_reader,
			reference_seq,
			"1", // Chromosome name
			0 // Minimum bridge length
		);
		
		// Generate the graph.
		std::vector <std::string> field_names_for_filter_by_assigned;
		generator.generate_graph(field_names_for_filter_by_assigned, true);
		
		// Check the graph.
		auto const &graph(generator.variant_graph());
		cmp.check_graph(graph, reference_seq);
	}
}


SCENARIO("Variant graph can be created correctly from a set of miscellaneous variants")
{
	GIVEN("A VCF file (1a)")
	{
		node_comparator cmp{
			{
				{0,		0,	0,	false,	"AAAA",	{}},
				{1,		4,	4,	true,	"A",	{{2, "G"}}},
				{2,		5,	5,	false,	"A",	{}},
				{3,		6,	6,	true,	"A",	{{4, "T"}, {4, "CC"}}},
				{4,		7,	8,	true,	"A",	{{5, "T"}, {5, "GGGG"}}},
				{5,		8,	12,	true,	"A",	{{7, "T"}}},
				{6,		9,	13,	false,	"A",	{{8, "CC"}}},
				{7,		10,	14,	false,	"A",	{{9, "GG"}}},
				{8,		11,	15,	false,	"A",	{}},
				{9,		12,	16,	false,	"AA",	{}},
				{10,	14,	18,	false,	"",		{}}
			}
		};
		
		test_variant_graph("test-1a.vcf", "test-1.fa", cmp);
	}
	
	GIVEN("A VCF file (1b)")
	{
		node_comparator cmp{
			{
				{0,		0,	0,	false,	"AAAA",	{}},
				{1,		4,	4,	true,	"A",	{{2, "G"}}},
				{2,		5,	5,	false,	"A",	{}},
				{3,		6,	6,	true,	"A",	{{4, "T"}, {4, "CC"}}},
				{4,		7,	8,	true,	"A",	{{5, "T"}, {5, "GGGG"}}},
				{5,		8,	12,	true,	"A",	{{7, "T"}}},
				{6,		9,	13,	false,	"A",	{{8, "CC"}}},
				{7,		10,	14,	false,	"A",	{{9, "GG"}}},
				{8,		11,	15,	false,	"A",	{}},
				{9,		12,	16,	false,	"AA",	{}},
				{10,	14,	18,	false,	"",		{}}
			}
		};
		
		test_variant_graph("test-1b.vcf", "test-1.fa", cmp);
	}
	
	GIVEN("A VCF file (2)")
	{
		node_comparator cmp{{
			{0,		0,	0,	true,	"GC",	{{4, "TTTT"}}},
			{1,		2,	2,	false,	"AA",	{{2, "C"}}},
			{2,		4,	4,	false,	"C",	{{3, "GG"}}},
			{3,		5,	6,	false,	"C",	{}},
			{4,		6,	7,	false,	"",		{}}
		}};
		test_variant_graph("test-2.vcf", "test-2.fa", cmp);
	}
	
	GIVEN("A VCF file (3)")
	{
		node_comparator cmp{{
			{0,		0,	0,	true,	"T",	{{10, "T"}}},
			{1,		1,	1,	false,	"GC",	{{8, "C"}}},
			{2,		3,	3,	false,	"TG",	{{3, "CCCC"}}},
			{3,		5,	7,	false,	"G",	{}},
			{4,		6,	8,	false,	"G",	{{5, "T"}}},
			{5,		7,	9,	false,	"AG",	{}},
			{6,		9,	11,	false,	"GC",	{{10, "TTTT"}}},
			{7,		11,	13,	false,	"A",	{{8, "G"}, {9, "C"}}},
			{8,		12,	14,	false,	"A",	{}},
			{9,		13,	15,	false,	"CC",	{}},
			{10,	15,	17,	false,	"",		{}}
		}};
		test_variant_graph("test-3.vcf", "test-3.fa", cmp);
	}
	
	GIVEN("A VCF file (4)")
	{
		node_comparator cmp{{
			{0,		0,	0,	true,	"T",	{{10, "T"}}},
			{1,		1,	1,	false,	"GC",	{{8, "C"}}},
			{2,		3,	3,	false,	"TG",	{{3, "CCCC"}}},
			{3,		5,	7,	false,	"G",	{}},
			{4,		6,	8,	false,	"G",	{{5, "T"}}},
			{5,		7,	9,	false,	"AG",	{}},
			{6,		9,	11,	false,	"GC",	{{10, "TTTT"}}},
			{7,		11,	13,	false,	"A",	{{8, "G"}, {9, "C"}, {9, ""}}},
			{8,		12,	14,	false,	"A",	{}},
			{9,		13,	15,	false,	"CC",	{}},
			{10,	15,	17,	false,	"GGGG",	{}},
			{11,	19,	21,	false,	"",		{}},
		}};
		test_variant_graph("test-4.vcf", "test-4.fa", cmp);
	}
}
