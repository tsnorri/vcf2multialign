/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <numeric>
#include <vcf2multialign/variant_format.hh>
#include "utility.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace {
	
	inline void setup_segment(v2m::aligned_segment const &src, v2m::aligned_segment &dst, std::size_t const pos_diff, bool const is_match)
	{
		if (is_match)
			dst.type = v2m::segment_type::MATCH;
		else
			dst.type = v2m::segment_type::MISMATCH;
		
		dst.aligned_position = src.aligned_position + pos_diff;
		dst.ref.position = src.ref.position + pos_diff;
		dst.alt.position = src.alt.position + pos_diff;
	}
}



namespace vcf2multialign {
	
	std::pair <std::uint16_t, std::uint16_t> count_set_genotype_values(vcf::variant const &var, std::uint16_t const alt_idx)
	{
		typedef std::pair <std::uint16_t, std::uint16_t> return_type;
		
		auto const &samples(var.samples());
		libbio_assert(!samples.empty());
		auto const *gt_field(v2m::get_variant_format(var).gt);
		libbio_assert(gt_field);
		auto const &first_sample(samples.front());
		auto const &gt((*gt_field)(first_sample)); // vector of sample_genotype
		
		// If 0 == alt_idx, count all non-zero values. Otherwise, count the instances of the given value.
		if (0 == alt_idx)
		{
			auto const count(std::accumulate(gt.begin(), gt.end(), std::int32_t(0), [](auto const sum, auto const &sample_gt) -> std::int32_t {
				return (0 == sample_gt.alt ? sum : 1 + sum);
			}));
			return return_type(count, gt.size());
		}
		else
		{
			auto const count(std::accumulate(gt.begin(), gt.end(), std::int32_t(0), [alt_idx](auto const sum, auto const &sample_gt) -> std::int32_t {
				return (alt_idx == sample_gt.alt ? 1 + sum : sum);
			}));
			return return_type(count, gt.size());
		}
	}
	
	
	void split_mixed_segment(aligned_segment &src, aligned_segment_vector &dst)
	{
		auto const ref_size(src.ref.string.size());
		auto const alt_size(src.alt.string.size()); // OK b.c. segment_type is MIXED.
		auto const head_size(alt_size <= ref_size ? alt_size : ref_size - 1);
		libbio_assert_gt(ref_size, 0);
		
		// Split the strings into head and tail parts like so:
		// Ref	GATTACA				GCTT
		// Alt	GCTT				GATTACA
		//		    ^tail begins	   ^tail begins
		// Then create smaller segments.
		
		// Considering the case where the alt start with a gap (segment_type::MIXED_ALT_STARTS_WITH_GAP == src.type),
		// the alt character position needs to be adjusted by one b.c. the parser sets the segment position to that of the first
		// alt character (which is set to the position of the previous non-gap character). Obviously this need not be done if
		// there actually are no alt characters.
		if (segment_type::MIXED_ALT_STARTS_WITH_GAP == src.type && alt_size)
			++src.alt.position;
		
		// Handle the head part.
		aligned_segment seg;
		{
			auto const head_range(
				rsv::zip(rsv::iota(std::size_t(0)), src.ref.string, src.alt.string)
				| rsv::take_exactly(head_size)
			);
			segment(
				head_range,
				[](auto const &tup) -> bool {								// Project
					auto const [i, refc, altc] = tup;
					return (refc == altc);
				},
				[&src, &seg](auto const &tup, auto const is_match){			// Start a segment
					auto const [i, refc, altc] = tup;
					setup_segment(src, seg, i, is_match);
				},
				[&src, &seg](auto const &tup, auto const is_match){			// Handle an item
					auto const [i, refc, altc] = tup;
					seg.ref.string += refc;
					if (!is_match)
						seg.alt.string += altc;
				},
				[&dst, &seg](auto const is_match){							// Finish a segment
					auto &new_seg(dst.emplace_back());
					check_segment(seg);
					// Move the strings.
					using std::swap;
					swap(new_seg, seg);
				}
			);
		}
		
		// Handle the tail part.
		{
			auto &ref_str(src.ref.string);
			auto &alt_str(src.alt.string);
			
			if (alt_size < ref_size)
			{
				ref_str.erase(ref_str.begin(), ref_str.begin() + alt_size);
				alt_str.clear();
			
				src.type = segment_type::DELETION;
				src.aligned_position += alt_size;
				src.ref.position += alt_size;
				src.alt.position += alt_size;
				check_segment(src);
			
				auto &new_seg(dst.emplace_back());
				using std::swap;
				swap(new_seg, src);
			}
			else if (ref_size < alt_size)
			{
				ref_str.erase(ref_str.begin(), ref_str.begin() + ref_size - 1);
				alt_str.erase(alt_str.begin(), alt_str.begin() + ref_size - 1);
			
				if (ref_str[0] == alt_str[0])
					src.type = segment_type::INSERTION;
				else
					src.type = segment_type::INSERTION_WITH_SNP;
			
				src.aligned_position += ref_size - 1;
				src.ref.position += ref_size - 1;
				src.alt.position += ref_size - 1;
				check_segment(src);
			
				auto &new_seg(dst.emplace_back());
				using std::swap;
				swap(new_seg, src);
			}
		}
	}
}
