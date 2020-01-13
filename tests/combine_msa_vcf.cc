/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <boost/format.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <catch2/catch.hpp>
#include <iterator>
#include <libbio/mmap_handle.hh>
#include "../combine-msa-vcf/combine_msa.hh"

namespace gen	= Catch::Generators;
namespace io	= boost::iostreams;
namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	inline std::string_view make_substring(std::string_view const &sv, std::string_view::const_iterator const &it)
	{
		auto const dist(std::distance(sv.begin(), it));
		return sv.substr(dist);
	}
	
	
	void compare_against_expected(std::string_view const &actual_sv, std::string_view const &expected_sv)
	{
		auto const res(std::mismatch(expected_sv.begin(), expected_sv.end(), actual_sv.begin(), actual_sv.end()));
		if (expected_sv.end() == res.first && actual_sv.end() == res.second)
			SUCCEED();
		else
		{
			INFO("First mismatch:");
			INFO("(expected): '" << make_substring(expected_sv, res.first)	<< '\'');
			INFO("(actual):   '" << make_substring(actual_sv, res.second)	<< '\'');
			FAIL();
		}
	}
	
	
	void compare_against_expected_path(std::string_view const &actual_sv, char const *expected_path)
	{
		lb::mmap_handle <char> expected_handle;
		expected_handle.open(expected_path);
		auto const &expected_sv(expected_handle.to_string_view());
		compare_against_expected(actual_sv, expected_sv);
	}
	
	
	void test_msa_merge(char const *test_data_dir, char const *expected_outcome, char const *ref_fa_name = "ref.fa", char const *alt_fa_name = "alt.fa", char const *expected_vcf_name = "expected.vcf")
	{
		auto fmt(boost::format("test-files/combine-msa-vcf/%s/%s"));
		auto const ref_path((fmt % test_data_dir % ref_fa_name).str());
		auto const alt_path((fmt % test_data_dir % alt_fa_name).str());
		
		WHEN("the files are merged")
		{
			std::string dst;
			io::filtering_ostream out(std::back_inserter(dst));
			v2m::combine_msa(ref_path.c_str(), alt_path.c_str(), nullptr, "chr1", 2, out);
			out.flush();
			
			THEN(expected_outcome)
			{
				auto const expected_path((fmt % test_data_dir % expected_vcf_name).str());
				compare_against_expected_path(dst, expected_path.c_str());
			}
		}
	}
}


SCENARIO("MSA combiner can handle empty sequences")
{
	GIVEN("An empty MSA")
	{
		test_msa_merge("empty", "the resulting VCF will have no variants", "empty.txt", "empty.txt", "empty.txt");
	}
}


SCENARIO("MSA combiner can merge sequences")
{
	GIVEN("A MSA with a SNP and a MNP")
	{
		test_msa_merge("msa-mnps", "the resulting VCF will have the variants");
	}
	
	GIVEN("A MSA with an insertion")
	{
		test_msa_merge("msa-ins", "the resulting VCF will have the insertion");
	}
	
	GIVEN("A MSA with a deletion")
	{
		test_msa_merge("msa-del", "the resulting VCF will have the deletion");
	}
	
	GIVEN("A MSA with an insertion and a SNP")
	{
		test_msa_merge("msa-ins-snp", "the resulting VCF will have the insertion and the SNP with the correct co-ordinates");
	}
	
	GIVEN("A MSA with insertions and deletions")
	{
		test_msa_merge("msa-indels-vcf-insertions-snps", "the resulting VCF will have the variants");
	}
	
	GIVEN("A MSA with insertions and deletions without description lines")
	{
		test_msa_merge("msa-indels-vcf-insertions-snps", "the resulting VCF will have the variants", "ref-bare.fa", "alt-bare.fa");
	}
}


SCENARIO("MSA combiner can merge sequences with gaps at start")
{
	GIVEN("A MSA with gaps at start (deletion)")
	{
		test_msa_merge("msa-gaps-at-start-del", "the resulting VCF will have the deletion in the first co-ordinate");
	}
	
	GIVEN("A MSA with gaps at start (insertion)")
	{
		test_msa_merge("msa-gaps-at-start-ins", "the resulting VCF will have the insertion (with a non-SNP) in the first co-ordinate");
	}
}


SCENARIO("MSA combiner can merge sequenes with mixed-type segments")
{
	GIVEN("A MSA with a mixed segment that reduces to a matching segment")
	{
		test_msa_merge("msa-mixed-1-1", "the resulting VCF will have no variants");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to non-matching segments")
	{
		test_msa_merge("msa-mixed-1-2", "the resulting VCF will have MNPs");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion without a SNP)")
	{
		test_msa_merge("msa-mixed-1-3", "the resulting VCF will have the expected variants");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion with a SNP)")
	{
		test_msa_merge("msa-mixed-1-4", "the resulting VCF will have the expected variants");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (deletion)")
	{
		test_msa_merge("msa-mixed-1-5", "the resulting VCF will have the expected variants");
	}
}


SCENARIO("MSA combiner can merge sequences and variants")
{
	GIVEN("A MSA with insertions and deletions and a VCF with insertions and SNPs")
	{
		char const *ref_path("test-files/combine-msa-vcf/msa-indels-vcf-insertions-snps/ref.fa");
		char const *alt_path("test-files/combine-msa-vcf/msa-indels-vcf-insertions-snps/alt.fa");
		char const *vcf_path("test-files/combine-msa-vcf/msa-indels-vcf-insertions-snps/vars.vcf");
		
		WHEN("the files are merged")
		{
			std::string dst;
			io::filtering_ostream out(std::back_inserter(dst));
			v2m::combine_msa(ref_path, alt_path, vcf_path, "chr1", 2, out);
			out.flush();
			
			THEN("the generated VCF matches the expected")
			{
				compare_against_expected_path(dst, "test-files/combine-msa-vcf/msa-indels-vcf-insertions-snps/expected-2.vcf");
			}
		}
	}
}
