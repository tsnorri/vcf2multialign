/*
 * Copyright (c) 2020-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <catch2/catch.hpp>
#include <libbio/fasta_reader.hh>
#include <vcf2multialign/variant_graph.hh>

namespace gen	= Catch::Generators;
namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	
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
	
	
	struct node_description
	{
		std::vector <alt_edge>	alt_edges;
		std::string				ref;
		std::size_t				node{};
		std::size_t				pos{};
		std::size_t				aln_pos{};
		
		node_description() = default;
		
		node_description(
			std::size_t node_,
			std::size_t pos_,
			std::size_t aln_pos_,
			std::string ref_,
			std::vector <alt_edge> alt_edges_
		):
			alt_edges(std::move(alt_edges_)),
			ref(std::move(ref_)),
			node(node_),
			pos(pos_),
			aln_pos(aln_pos_)
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
		
		void check_graph(v2m::sequence_type const &reference, v2m::variant_graph const &graph)
		{
			v2m::variant_graph_walker walker(reference, graph);
			for (auto const &desc : m_node_descriptions)
			{
				REQUIRE(walker.advance());
				CHECK(desc.node == walker.node());
				INFO("Node " << desc.node)
				CHECK(desc.pos == walker.ref_position());
				INFO("Ref pos " << desc.pos);
				CHECK(desc.aln_pos == walker.aligned_position());
				INFO("Aln pos " << desc.aln_pos);
				
				auto const &expected_alt_edges(desc.alt_edges);
				auto const &actual_alt_edges(walker.alt_edges());
				CHECK(expected_alt_edges.size() == actual_alt_edges.size());
				for (auto const &[edge_idx, tup] : rsv::enumerate(rsv::zip(expected_alt_edges, actual_alt_edges)))
				{
					auto const &[expected_edge, actual_edge] = tup;
					auto const &[actual_target, actual_label] = actual_edge;
					INFO("ALT edge " << edge_idx);
					INFO("Expected target: " << expected_edge.target_node);
					INFO("Actual target:   " << actual_target);
					INFO("Expected label:  " << expected_edge.label);
					INFO("Actual label:    " << actual_label);
					CHECK(expected_edge.target_node == actual_target);
					CHECK(expected_edge.label == actual_label);
				}
			}
			
			REQUIRE(!walker.advance());
		}
	};
	
	
	void test_variant_graph(char const *vcf_name, char const *fasta_name, node_comparator &cmp)
	{
		std::string const data_dir("test-files/variant-graph/");
		auto const vcf_path(data_dir + vcf_name);
		auto const fasta_path(data_dir + fasta_name);
		
		v2m::sequence_type ref_seq;
		REQUIRE(lb::read_single_fasta_sequence(fasta_path.c_str(), ref_seq));
		
		v2m::variant_graph graph;
		v2m::build_graph_statistics stats;
		v2m::build_variant_graph(ref_seq, vcf_path.c_str(), "1", graph, stats, nullptr);
		
		cmp.check_graph(ref_seq, graph);
	}
}


SCENARIO("Variant graph can be created correctly from a set of miscellaneous variants")
{
	GIVEN("A VCF file (1a)")
	{
		node_comparator cmp{
			{
				{0,		0,	0,	"AAAA",	{}},
				{1,		4,	4,	"A",	{{2, "G"}}},
				{2,		5,	5,	"A",	{}},
				{3,		6,	6,	"A",	{{4, "T"}, {4, "CC"}}},
				{4,		7,	8,	"A",	{{5, "T"}, {5, "GGGG"}}},
				{5,		8,	12,	"A",	{{7, "T"}}},
				{6,		9,	13,	"A",	{{8, "CC"}}},
				{7,		10,	14,	"A",	{{9, "GG"}}},
				{8,		11,	15,	"A",	{}},
				{9,		12,	16,	"AA",	{}},
				{10,	14,	18,	"",		{}}
			}
		};
		
		test_variant_graph("test-1a.vcf", "test-1.fa", cmp);
	}
	
	GIVEN("A VCF file (1b)")
	{
		node_comparator cmp{
			{
				{0,		0,	0,	"AAAA",	{}},
				{1,		4,	4,	"A",	{{2, "G"}}},
				{2,		5,	5,	"A",	{}},
				{3,		6,	6,	"A",	{{4, "T"}, {4, "CC"}}},
				{4,		7,	8,	"A",	{{5, "T"}, {5, "GGGG"}}},
				{5,		8,	12,	"A",	{{7, "T"}}},
				{6,		9,	13,	"A",	{{8, "CC"}}},
				{7,		10,	14,	"A",	{{9, "GG"}}},
				{8,		11,	15,	"A",	{}},
				{9,		12,	16,	"AA",	{}},
				{10,	14,	18,	"",		{}}
			}
		};
		
		test_variant_graph("test-1b.vcf", "test-1.fa", cmp);
	}
	
	GIVEN("A VCF file (2)")
	{
		node_comparator cmp{{
			{0,		0,	0,	"GC",	{{4, "TTTT"}}},
			{1,		2,	2,	"AA",	{{2, "C"}}},
			{2,		4,	4,	"C",	{{3, "GG"}}},
			{3,		5,	6,	"C",	{}},
			{4,		6,	7,	"",		{}}
		}};
		test_variant_graph("test-2.vcf", "test-2.fa", cmp);
	}
	
	GIVEN("A VCF file (3)")
	{
		node_comparator cmp{{
			{0,		0,	0,	"T",	{{10, "T"}}},
			{1,		1,	1,	"GC",	{{8, "C"}}},
			{2,		3,	3,	"TG",	{{3, "CCCC"}}},
			{3,		5,	7,	"G",	{}},
			{4,		6,	8,	"G",	{{5, "T"}}},
			{5,		7,	9,	"AG",	{}},
			{6,		9,	11,	"GC",	{{10, "TTTT"}}},
			{7,		11,	13,	"A",	{{8, "G"}, {9, "C"}}},
			{8,		12,	14,	"A",	{}},
			{9,		13,	15,	"CC",	{}},
			{10,	15,	17,	"",		{}}
		}};
		test_variant_graph("test-3.vcf", "test-3.fa", cmp);
	}
	
	GIVEN("A VCF file (4)")
	{
		node_comparator cmp{{
			{0,		0,	0,	"T",	{{10, "T"}}},
			{1,		1,	1,	"GC",	{{8, "C"}}},
			{2,		3,	3,	"TG",	{{3, "CCCC"}}},
			{3,		5,	7,	"G",	{}},
			{4,		6,	8,	"G",	{{5, "T"}}},
			{5,		7,	9,	"AG",	{}},
			{6,		9,	11,	"GC",	{{10, "TTTT"}}},
			{7,		11,	13,	"A",	{{8, "G"}, {9, "C"}, {9, ""}}},
			{8,		12,	14,	"A",	{}},
			{9,		13,	15,	"CC",	{}},
			{10,	15,	17,	"GGGG",	{}},
			{11,	19,	21,	"",		{}},
		}};
		test_variant_graph("test-4.vcf", "test-4.fa", cmp);
	}
}
