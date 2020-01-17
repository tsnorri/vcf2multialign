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
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
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
	
	
	void combine_msa_and_compare_against_expected_path(
		v2m::vector_type const &ref_seq,
		v2m::vector_type const &alt_seq,
		char const *vcf_path,
		char const *expected_outcome,
		std::string const &expected_path
	)
	{
		std::string dst;
		io::filtering_ostream out(std::back_inserter(dst));
		v2m::combine_msa(ref_seq, alt_seq, vcf_path, "chr1", 2, out, false);
		out.flush();
			
		THEN(expected_outcome)
		{
			compare_against_expected_path(dst, expected_path.c_str());
		}
	}
	
	
	void test_msa_vcf_merge(char const *test_data_dir, char const *expected_outcome, char const *ref_fa_name = "ref.fa", char const *alt_fa_name = "alt.fa", char const *vcf_name = "vars.vcf", char const *expected_vcf_name = "expected.vcf")
	{
		// Read the input.
		v2m::vector_type ref_seq, alt_seq;
		auto fmt(boost::format("test-files/combine-msa-vcf/%s/%s"));
		auto const ref_path((fmt % test_data_dir % ref_fa_name).str());
		auto const alt_path((fmt % test_data_dir % alt_fa_name).str());
		auto const vcf_path(vcf_name ? (fmt % test_data_dir % vcf_name).str() : "");
		auto const vcf_path_c(vcf_name ? vcf_path.c_str() : nullptr);
		auto const expected_path((fmt % test_data_dir % expected_vcf_name).str());
		v2m::read_single_fasta_seq(ref_path.c_str(), ref_seq, nullptr, false);
		v2m::read_single_fasta_seq(alt_path.c_str(), alt_seq, nullptr, false);
		
		INFO("Ref path:          " << ref_path);
		INFO("Alt path:          " << alt_path);
		INFO("VCF path:          " << vcf_path);
		INFO("Expected VCF path: " << expected_path);

		WHEN("the files are merged")
		{
			combine_msa_and_compare_against_expected_path(ref_seq, alt_seq, vcf_path_c, expected_outcome, expected_path);
		}
		
		// Repeat the test with the sequences converted to lowercase.
		std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), [](auto const c){
			return std::tolower(c);
		});
		std::transform(alt_seq.begin(), alt_seq.end(), alt_seq.begin(), [](auto const c){
			return std::tolower(c);
		});
		
		WHEN("the files are merged with lowercase FASTA inputs")
		{
			combine_msa_and_compare_against_expected_path(ref_seq, alt_seq, vcf_path_c, expected_outcome, expected_path);
		}
	}
	
	
	void test_msa_merge(char const *test_data_dir, char const *expected_outcome, char const *ref_fa_name = "ref.fa", char const *alt_fa_name = "alt.fa", char const *expected_vcf_name = "expected.vcf")
	{
		test_msa_vcf_merge(test_data_dir, expected_outcome, ref_fa_name, alt_fa_name, nullptr, expected_vcf_name);
	}
	
	
	template <typename t_item>
	inline std::size_t distance_to_begin(std::vector <t_item> const &vec, typename std::vector <t_item>::const_iterator it)
	{
		return std::distance(vec.begin(), it);
	}
}


SCENARIO("MSA combiner can locate the range of overlapping segments")
{
	struct variant_mock
	{
		std::size_t	position{};
		
		variant_mock() = default;
		variant_mock(std::size_t position_): position(position_) {}
		std::size_t zero_based_pos() const { return position; }
	};
	
	struct variant_record_mock
	{
		variant_mock	variant;
		std::size_t		size{};
		
		variant_record_mock() = default;
		variant_record_mock(std::size_t const position, std::size_t const size_):
			variant(position),
			size(size_)
		{
		}
	};
	
	GIVEN("A vector of segments with two deletions")
	{
		// Corresponds to:
		// GATTACAGATTACA
		// GA---CA------A
		
		std::vector <v2m::aligned_segment> const vec({
			{"GA",		"",	0,	0,	0,	v2m::segment_type::MATCH},
			{"TTA",		"",	2,	1,	2,	v2m::segment_type::DELETION},
			{"CA",		"", 5,	2,	5,	v2m::segment_type::MATCH},
			{"GATTAC",	"", 7,	3,	7,	v2m::segment_type::DELETION},
			{"A",		"", 13,	4,	13,	v2m::segment_type::MATCH}
		});
		
		WHEN("a variant that overlaps with the third segment at the boundary of a preceding deletion is tested")
		{
			variant_record_mock const var(2, 1);
			
			THEN("a range with the third segment is returned")
			{
				auto const it_pair(find_overlapping_segment_range(vec.begin(), vec.end(), var));
				CHECK(2 == distance_to_begin(vec, it_pair.first));
				CHECK(3 == distance_to_begin(vec, it_pair.second));
			}
		}
		
		WHEN("a variant that overlaps with the third segment at the boundary of a succeeding deletion is tested")
		{
			variant_record_mock const var(3, 1);
			
			THEN("a range with the third and the fourth segment is returned")
			{
				auto const it_pair(find_overlapping_segment_range(vec.begin(), vec.end(), var));
				CHECK(2 == distance_to_begin(vec, it_pair.first));
				// Since the deletion occurs right after the overlapping segment, the variant “overlaps” with the
				// gap characters of the deletion.
				CHECK(4 == distance_to_begin(vec, it_pair.second));
			}
		}
		
		WHEN("a variant that overlaps with the first three segments is tested")
		{
			variant_record_mock var(1, 2);
			
			THEN("a range with the first three segments is returned")
			{
				auto const it_pair(find_overlapping_segment_range(vec.begin(), vec.end(), var));
				CHECK(0 == distance_to_begin(vec, it_pair.first));
				CHECK(3 == distance_to_begin(vec, it_pair.second));
			}
		}
		
		WHEN("a variant that overlaps with all of the segments (except the first character of the first segement) is tested")
		{
			variant_record_mock var(1, 4);
			
			THEN("a range with all of the segments is returned")
			{
				auto const it_pair(find_overlapping_segment_range(vec.begin(), vec.end(), var));
				CHECK(0 == distance_to_begin(vec, it_pair.first));
				CHECK(5 == distance_to_begin(vec, it_pair.second));
			}
		}
	}
}


SCENARIO("MSA combiner can handle empty sequences")
{
	GIVEN("An empty MSA")
	{
		test_msa_merge("empty", "the resulting VCF will have no variants", "empty.txt", "empty.txt", "expected.vcf");
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
	
	GIVEN("A MSA with an insertion")
	{
		test_msa_merge("msa-ins-vcf-multiple-overlaps", "the resulting VCF will have the insertions");
	}
	
	GIVEN("A MSA with a deletion")
	{
		test_msa_merge("msa-del", "the resulting VCF will have the deletion");
	}
	
	GIVEN("A MSA with two deletions")
	{
		test_msa_merge("msa-del-vcf-overlaps", "the resulting VCF will have the deletions");
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
	
	GIVEN("A MSA with different types of segments (to test the remaining transitions)")
	{
		test_msa_merge("msa-different-segments", "the resulting VCF will have the variants");
	}
	
	GIVEN("A MSA with a mixed segment in the beginning (gap in both)")
	{
		test_msa_merge("msa-mixed-beginning-1", "the resulting VCF will have the variants");
	}
	
	GIVEN("A MSA with a mixed segment in the beginning (gap in ref)")
	{
		test_msa_merge("msa-mixed-beginning-2", "the resulting VCF will have the variants");
	}
	
	GIVEN("A MSA with a mixed segment in the beginning (gap in alt)")
	{
		test_msa_merge("msa-mixed-beginning-3", "the resulting VCF will have the variants");
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


SCENARIO("MSA combiner can merge sequences with mixed-type segments (ref precedes)")
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


SCENARIO("MSA combiner can merge sequences with mixed-type segments (ref precedes, less gaps)")
{
	GIVEN("A MSA with a mixed segment that reduces to a matching segment")
	{
		test_msa_merge("msa-mixed-1-1", "the resulting VCF will have no variants", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to non-matching segments")
	{
		test_msa_merge("msa-mixed-1-2", "the resulting VCF will have MNPs", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion without a SNP)")
	{
		test_msa_merge("msa-mixed-1-3", "the resulting VCF will have the expected variants", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion with a SNP)")
	{
		test_msa_merge("msa-mixed-1-4", "the resulting VCF will have the expected variants", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (deletion)")
	{
		test_msa_merge("msa-mixed-1-5", "the resulting VCF will have the expected variants", "ref-2.fa", "alt-2.fa");
	}
}


SCENARIO("MSA combiner can merge sequences with mixed-type segments (alt precedes)")
{
	GIVEN("A MSA with a mixed segment that reduces to a matching segment in the middle")
	{
		test_msa_merge("msa-mixed-2-1", "the resulting VCF will have a MNP and a deletion");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to non-matching segments")
	{
		test_msa_merge("msa-mixed-2-2", "the resulting VCF will have MNPs and a deletion");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion without a SNP)")
	{
		test_msa_merge("msa-mixed-2-3", "the resulting VCF will have the expected variants");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion with a SNP)")
	{
		test_msa_merge("msa-mixed-2-4", "the resulting VCF will have the expected variants");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (deletion)")
	{
		test_msa_merge("msa-mixed-2-5", "the resulting VCF will have the expected variants");
	}
}


SCENARIO("MSA combiner can merge sequences with mixed-type segments (alt precedes, less gaps)")
{
	GIVEN("A MSA with a mixed segment that reduces to a matching segment in the middle")
	{
		test_msa_merge("msa-mixed-2-1", "the resulting VCF will have a MNP and a deletion", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to non-matching segments")
	{
		test_msa_merge("msa-mixed-2-2", "the resulting VCF will have MNPs and a deletion", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion without a SNP)")
	{
		test_msa_merge("msa-mixed-2-3", "the resulting VCF will have the expected variants", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (insertion with a SNP)")
	{
		test_msa_merge("msa-mixed-2-4", "the resulting VCF will have the expected variants", "ref-2.fa", "alt-2.fa");
	}
	
	GIVEN("A MSA with a mixed segment that reduces to different types of segments (deletion)")
	{
		test_msa_merge("msa-mixed-2-5", "the resulting VCF will have the expected variants", "ref-2.fa", "alt-2.fa");
	}
}


SCENARIO("MSA combiner can merge sequences with mixed-type segemnts (miscellaneous)")
{
	GIVEN("A MSA with a mixed segment that reduces to a matching segment, a mismatching segment and a matching segment")
	{
		test_msa_merge("msa-mixed-multiple-types", "the resulting VCF will have the expected variants");
	}
}


SCENARIO("MSA combiner can merge sequences and variants")
{
	GIVEN("A MSA with insertions and deletions and a VCF with insertions and SNPs")
	{
		test_msa_vcf_merge("msa-indels-vcf-insertions-snps", "the resulting VCF will have the expected variants", "ref.fa", "alt.fa", "vars.vcf", "expected-2.vcf");
	}
	
	GIVEN("An identity MSA and a VCF with overlapping variants")
	{
		test_msa_vcf_merge("msa-snps-vcf-overlaps", "the resulting VCF will have the expected variants", "ref.fa", "ref.fa");
	}
	
	GIVEN("A MSA with SNPs and a VCF with overlapping variants")
	{
		test_msa_vcf_merge("msa-snps-vcf-overlaps", "the resulting VCF will have the expected variants", "ref.fa", "alt.fa", "vars-2.vcf", "expected-2.vcf");
	}
	
	GIVEN("A MSA with two deletions and two variants that overlap with the deletions")
	{
		test_msa_vcf_merge("msa-del-vcf-overlaps", "the resulting VCF will have the deletions", "ref.fa", "alt.fa", "vars.vcf", "expected-2.vcf");
	}
	
	GIVEN("A MSA with an insertion and overlapping variants")
	{
		test_msa_vcf_merge("msa-ins-vcf-multiple-overlaps", "the resulting VCF will have the insertions", "ref.fa", "alt.fa", "vars.vcf", "expected-2.vcf");
	}
}
