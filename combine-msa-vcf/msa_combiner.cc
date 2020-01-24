/*
 * Copyright (c) 2019–2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/variant_printer.hh>
#include <range/v3/all.hpp>
#include <vcf2multialign/variant_format.hh>
#include "algorithms.hh"
#include "msa_combiner.hh"
#include "utility.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {

	// Determine the variant position relative to segment start.
	inline std::size_t segment_relative_variant_position(v2m::aligned_segment const &seg, v2m::variant_record const &var)
	{
		auto const var_pos(var.variant.zero_based_pos());
		auto const alt_pos(seg.alt.position);
		libbio_assert_lte(alt_pos, var_pos);
		return var_pos - alt_pos;
	}


	template <typename t_item>
	struct print_helper
	{
		std::vector <t_item> const	*vector{};

		print_helper(std::vector <t_item> const &vector_): vector(&vector_) {}
	};

	template <typename t_item>
	std::ostream &operator<<(std::ostream &os, print_helper <t_item> const &ph)
	{
		std::size_t i(0);
		for (auto const &item : *ph.vector)
		{
			os << "[" << i << "]: " << item << '\n';
			++i;
		}
		return os;
	}
}


namespace vcf2multialign {

	void msa_combiner::process_one_segment_msa(
		aligned_segment const &seg,
		std::int32_t overlap_count,
		overlap_counter::const_iterator overlap_it,
		overlap_counter::const_iterator const overlap_end
	)
	{
		// Output S/MNPs.
		switch (seg.type)
		{
			case segment_type::MATCH:
				return;
				
			case segment_type::MISMATCH:
			{
				std::string_view const ref_sv(seg.ref.string);
				std::string_view const alt_sv(seg.alt.string);

				auto const overlap_range(ranges::subrange(overlap_it, overlap_end));
				std::size_t begin_idx{};
				for (auto const &oc : overlap_range)
				{
					auto const end_idx(std::min(alt_sv.size(), std::size_t(oc.position - seg.alt.position)));
					auto const &desc(m_variant_writer.emplace_back(
						variant_description(
							seg.ref.position + begin_idx,
							ref_sv.substr(begin_idx, end_idx - begin_idx),
							alt_sv.substr(begin_idx, end_idx - begin_idx),
							m_ploidy,
							overlap_count,
							variant_origin::MSA
						)
					));
					libbio_assert_lt(0, desc.ref.size());
					
					begin_idx = end_idx;
					overlap_count = oc.running_sum;
				}
				
				if (begin_idx < alt_sv.size())
				{
					auto const &desc(m_variant_writer.emplace_back(
						variant_description(
							seg.ref.position + begin_idx,
							ref_sv.substr(begin_idx),
							alt_sv.substr(begin_idx),
							m_ploidy,
							overlap_count,
							variant_origin::MSA
						)
					));
					libbio_assert_lt(0, desc.ref.size());
				}
				break;
			}
				
			case segment_type::DELETION:
			{
				libbio_assert_eq(overlap_it, overlap_end);
				auto const &desc(m_variant_writer.emplace_back(
					variant_description(
						seg.ref.position,
						seg.ref.string,
						"",
						m_ploidy,
						overlap_count,
						variant_origin::MSA
					)
				));
				libbio_assert_lt(0, desc.ref.size());
				break;
			}
				
			case segment_type::INSERTION:
			case segment_type::INSERTION_WITH_SNP:
			{
				auto const overlap_range(ranges::subrange(overlap_it, overlap_end));
				// Determine the max. overlap count.
				auto const max_overlaps_in_range(
					ranges::empty(overlap_range)
					? std::int32_t(0)
					: ranges::max(overlap_range | rsv::transform([](auto const &oc){ return oc.running_sum; }))
				);
				auto const max_overlaps(std::max(overlap_count, max_overlaps_in_range));
				auto const &desc(m_variant_writer.emplace_back(
					variant_description(
						seg.ref.position,
						seg.ref.string,
						seg.alt.string,
						m_ploidy,
						max_overlaps,
						variant_origin::MSA
					)
				));
				libbio_assert_lt(0, desc.ref.size());
				break;
			}
				
			case segment_type::MIXED:
				throw std::runtime_error("Unexpected segment type");
		}
	}
	
	
	// Generate variants from MSA only.
	std::size_t msa_combiner::process_variants_msa(
		std::size_t const max_alt_pos,
		aligned_segment_vector::const_iterator seg_it,
		aligned_segment_vector::const_iterator const seg_end,
		class overlap_counter const &overlap_counter
	)
	{
		if (seg_it == seg_end) return 0;
		
		// Determine the overlap count range for the first segment.
		auto const &first_seg(*seg_it);
		auto const first_seg_alt_pos(first_seg.alt.position);
		auto const first_seg_alt_end_pos(first_seg.alt_end());
		auto overlap_it(overlap_counter.find_initial(first_seg_alt_pos));
		
		// Add the variants from the MSA to m_variant_writer.
		// This will result in two partitions of sorted variants.
		std::size_t handled_count(0);
		while (seg_it != seg_end)
		{
			// Check that max_alt_pos does not overlap with the current segment.
			auto const &seg(*seg_it);
			auto const seg_alt_pos(seg.alt.position);
			auto const seg_alt_end_pos(seg.alt_end());
			if (max_alt_pos < seg_alt_end_pos)
				break;
			
			overlap_counter.update_overlap_iterator_if_needed(overlap_it, seg_alt_pos);
			auto const initial_overlap_count(overlap_counter.initial_count(overlap_it));
			auto const overlap_end(overlap_counter.find_end(overlap_it, seg_alt_end_pos));
			process_one_segment_msa(seg, initial_overlap_count, overlap_it, overlap_end);
			
			overlap_it = overlap_end;
			++handled_count;
			++seg_it;
		}
		
		return handled_count;
	}
	
	
	// Generate variants from the given variant_record associated with overlapping aligned_segments.
	void msa_combiner::process_variant_in_range(
		variant_record const &var,
		aligned_segment_vector::const_iterator aligned_segment_begin,
		aligned_segment_vector::const_iterator const aligned_segment_end,
		std::int32_t const max_overlaps
	)
	{
		// The overlap types may be classified as follows.
		// (1) and (2) are mutually exclusive. If a variant is split, it could be marked as phased
		// with a phase set identifier. Aligned_segments have already been classified to 
		// MATCH, MISMATCH, DELETION, INSERTION and INSERTION_WITH_SNP.
		//
		// Head / tail
		// ===========
		// 				Seg_1 …	Seg_n
		// Ref / Alt	XXXXX	XXXXX
		// Var			   YY   YY
		//				^^^		  ^^^
		// These positions may be handled by considering the MSA only.
		//
		// 1. Variant overlaps with alt polymorphism
		// =========================================
		//
		// Ref	GATTACA		GAT
		// Alt	TCAACGC		GGGTACA
		// Var	  GG		 CC
		//
		// In this case, only the REF in the variant needs to be rewritten.
		// The variant could be segmented, though, by removing the parts that
		// match the original reference.
		// The GT field could be taken directly from the variant.
		//
		// 2. Variant overlaps with alt insertion
		// ======================================
		//
		// Ref	GAA			GAT
		// Alt	GATTACA		GATTACA
		// Var	    C		   CC
		// Var2	  TTC		  TCC	<- Variants need to be rewritten like this s.t. they overlap with REF.
		//
		// Variant overlaps with deletion
		// ==============================
		//
		// Ref	GATTACA
		// Alt	G-----G
		// Var	C-----C
		//
		// The REF field in the variant needs to be rewritten, in this case from GG to GATTACA.
		//
		// Multiple ALTs
		// =============
		// In case a variant contains multiple ALTs, handle each ALT separately, push the associated
		// variant information to a stack and after handling either the position in question or the
		// set of variants that overlap with each other, filter or rewrite the stack contents.
		
		libbio_assert_neq(aligned_segment_begin, aligned_segment_end);
		libbio_assert_lt(0, max_overlaps);
		
		auto const &first_seg(*aligned_segment_begin);
		libbio_assert_neq(first_seg.type, segment_type::DELETION);
		auto const var_pos_seg_relative(segment_relative_variant_position(first_seg, var));
		
		// Handle the variant.
		{
			auto const range(ranges::subrange(aligned_segment_begin, aligned_segment_end));
			std::size_t alt_idx{};
			for (auto const &var_alt : var.variant.alts())
			{
				++alt_idx;
				// Skip if the ALT allele was unknown.
				if ("*" == var_alt.alt)
					continue;

				auto const [gt_count, ploidy] = count_set_genotype_values(var.variant, alt_idx);
				libbio_always_assert_eq_msg(ploidy, m_ploidy, "Line ", var.variant.lineno(), ": expected the sample ploidy to match the passed value, got ", ploidy, '.');
				// Skip if no GT values were set.
				if (0 == gt_count)
					continue;

				auto left_pad(var_pos_seg_relative);
				std::size_t var_ref_characters_remaining(var.size); // Number of remaining variant reference characters, i.e. in seg.alt.
				std::size_t var_alt_characters_remaining(var_alt.alt.size());
				std::size_t total_alt_characters_consumed{};

				auto &desc(m_variant_writer.emplace_back(variant_description(m_ploidy, gt_count, variant_origin::VC)));
				desc.overlap_count = max_overlaps - 1;
				std::string_view const var_alt_sv(var_alt.alt);
				
				switch (first_seg.type)
				{
					case segment_type::MISMATCH:
					case segment_type::MATCH:
						desc.position = first_seg.ref.position + var_pos_seg_relative;
						break;
						
					case segment_type::INSERTION:
					case segment_type::INSERTION_WITH_SNP:
						desc.position = first_seg.ref.position;
						break;
						
					case segment_type::MIXED:
					case segment_type::DELETION:
						throw std::runtime_error("Unexpected segment type");
				}
				
				for (auto const &seg : range)
				{
					switch (seg.type)
					{
						case segment_type::MISMATCH:
						case segment_type::MATCH:
						{
							libbio_assert_lt(0, var_ref_characters_remaining);
							
							std::string_view const ref_sv(seg.ref.string);
							std::string_view const alt_sv(segment_type::MATCH == seg.type ? seg.ref.string : seg.alt.string);
							libbio_assert_eq(ref_sv.size(), alt_sv.size());
							
							auto const tail_length(alt_sv.size() - left_pad);
							auto const ref_characters_available(std::min(var_ref_characters_remaining, tail_length));
							auto const alt_characters_available(std::min(var_alt_characters_remaining, tail_length));
							auto const var_alt_start(var_alt_sv.size() - var_alt_characters_remaining); // Take into account the characters appended in the previous segments.
							desc.ref += ref_sv.substr(left_pad, ref_characters_available);
							desc.alt += var_alt_sv.substr(var_alt_start, alt_characters_available);
							
							var_ref_characters_remaining -= ref_characters_available;
							var_alt_characters_remaining -= alt_characters_available;
							break;
						}
						
						case segment_type::DELETION:
						{
							libbio_assert_eq(0, left_pad);
							
							// If we are at the end of the range, don’t append any characters as it would cause
							// an unnecessary overlap with the variant created from the deletion segment.
							// Continue the loop, though, in case assertions have been enabled.
							if (0 < var_ref_characters_remaining)
							{
								// Append the ref characters to the currently handled variant.
								desc.ref += seg.ref.string;
							}
							break;
						}
						
						case segment_type::INSERTION:
						case segment_type::INSERTION_WITH_SNP:
						{
							libbio_assert_lt(0, var_ref_characters_remaining);

							std::string_view const alt_sv(seg.alt.string);
							
							// Add the pad characters if needed.
							desc.alt += alt_sv.substr(0, left_pad);
							
							// Add the remaining characters.
							auto const tail_length(alt_sv.size() - left_pad);
							auto const ref_characters_available(std::min(var_ref_characters_remaining, tail_length));
							auto const alt_characters_available(std::min(var_alt_characters_remaining, tail_length));
							auto const var_alt_start(var_alt_sv.size() - var_alt_characters_remaining);
							desc.ref += seg.ref.string; // var_ref_characters_remaining is at least one.
							desc.alt += var_alt_sv.substr(var_alt_start, alt_characters_available);
							
							var_ref_characters_remaining -= ref_characters_available;
							var_alt_characters_remaining -= alt_characters_available;
							break;
						}
						
						case segment_type::MIXED:
							throw std::runtime_error("Unexpected segment type");
					}
				
					left_pad = 0;
				}
				
				// Append the remaining characters if needed.
				libbio_assert_eq(0, var_ref_characters_remaining);
				desc.alt += var_alt_sv.substr(var_alt_sv.size() - var_alt_characters_remaining);

				libbio_assert_lt(0, desc.ref.size());
			}
		}
	}
	
	
	// Handle a number of segments with variants that overlap with them.
	void msa_combiner::process_variants(msa_segmentation &segmentation)
	{
		// Consider a situation like this:
		//
		// Ref	G
		// Alt	GATTACA
		// Var	  G
		// Var	    T
		// Var	  GCC
		//
		// Suppose G, T and GCC are ALT column values. Since Alt (the ad-hoc reference) by itself
		// represents one chromosome copy, G -> GATTACA is one variant. Considering a diploid individual,
		// either G, T, G and T or GCC could be applied to it to determine the variant in the other copy.
		// In polyploid individuals, even more variants could be applied.
		//
		// The current implementation converts all the variants reported by the VC to the reference co-ordinates
		// and creates variant records from the MSA while taking into account the boundaries of the variants
		// reported by the VC. The number of overlapping variants affects the GT values that are set to them.
		
		segmentation.overlap_counter.update_running_sums();
		auto const initial_variant_count(m_variant_writer.size());
		if (segmentation.overlapping_variants.empty())
		{
			// Handle MSA only.
			auto const end(segmentation.overlapping_segments.end());
			auto const handled_count(
				process_variants_msa(
					SIZE_MAX,
					segmentation.overlapping_segments.begin(),
					end,
					segmentation.overlap_counter
				)
			);
			segmentation.overlapping_segments.erase(
				segmentation.overlapping_segments.begin(),
				segmentation.overlapping_segments.begin() + handled_count
			);
			
			m_variant_writer.merge_output_variants(initial_variant_count);
			m_variant_writer.filter_processed_variants_and_output(SIZE_MAX);
		}
		else
		{
			libbio_assert(!segmentation.overlapping_segments.empty());
			libbio_assert_neq(SIZE_MAX, segmentation.overlapping_segments.front().alt.position);
			
			std::size_t min_unhandled_alt_pos(SIZE_MAX);
			std::size_t min_unhandled_ref_pos(SIZE_MAX);
			
			auto var_it(segmentation.overlapping_variants.begin());
			auto seg_it(segmentation.overlapping_segments.begin());
			auto overlap_it(segmentation.overlap_counter.begin());
			auto const var_end(segmentation.overlapping_variants.end());
			auto const seg_end(segmentation.overlapping_segments.end());
			auto const overlap_end(segmentation.overlap_counter.end());
			
			while (var_it != var_end)
			{
				auto const &var(*var_it);
				auto const var_pos(var.variant.zero_based_pos());
				auto const var_end_pos(var_pos + var.size);
				
				// Determine the segments that overlap with var.
				auto const seg_range(find_overlapping_segment_range(seg_it, seg_end, var));
				libbio_assert_neq_msg(seg_range.first, seg_end, "There should be at least one element not less than (overlapping with) the variant.");
				seg_it = seg_range.first;
				auto const seg_current_end(seg_range.second);
				
				// If the last segment satisfies the condition, the variant may extend over the end of the segment.
				if (seg_end == seg_current_end && (seg_current_end - 1)->alt_end() < var.variant.zero_based_pos() + var.size)
				{
					min_unhandled_ref_pos = std::min(min_unhandled_ref_pos, seg_it->ref.position); // Use the segment start b.c. converting the variant co-ordinate is difficult.
					min_unhandled_alt_pos = std::min(min_unhandled_alt_pos, var_pos);
					var_it->is_skipped = true;
					++var_it;
					continue;
				}
				
				// Determine the max. overlap count.
				auto const res(segmentation.overlap_counter.max_overlaps_in_range(overlap_it, overlap_end, var_pos, var_end_pos));
				auto const max_overlaps(res.first);
				libbio_assert_msg(
					0 < max_overlaps || 0 == count_set_genotype_values(var_it->variant, 0).first,
					"Got ", max_overlaps, " for maximum overlaps. Overlap list:\n", print_helper(segmentation.overlap_counter.overlap_counts()),
					"range: [",
					std::distance(segmentation.overlap_counter.overlap_counts().begin(), overlap_it), ", ",
					std::distance(segmentation.overlap_counter.overlap_counts().begin(), overlap_end), ") ",
					"var_pos: ", var_pos, " var_end_pos: ", var_end_pos
				);
				overlap_it = res.second;
				
				// Process the variant and the found segments.
				if (0 < max_overlaps)
					process_variant_in_range(var, seg_it, seg_current_end, max_overlaps);
				++var_it;
			}
			
			// Erase the handled variants.
			segmentation.overlapping_variants.erase(
				std::remove_if(
					segmentation.overlapping_variants.begin(),
					segmentation.overlapping_variants.end(),
					[](auto const &var){
						return !var.is_skipped;
					}
				),
				segmentation.overlapping_variants.end()
			);
			// Reset the flag.
			for (auto &var : segmentation.overlapping_variants)
				var.is_skipped = false;
			
			// Maintain the sorted order.
			m_variant_writer.merge_output_variants(initial_variant_count);
			
			auto const new_variant_count(m_variant_writer.size());
			auto const handled_count(
				process_variants_msa(
					min_unhandled_alt_pos,
					segmentation.overlapping_segments.begin(),
					segmentation.overlapping_segments.end(),
					segmentation.overlap_counter
				)
			);
			// Erase the handled segments.
			segmentation.overlapping_segments.erase(
				segmentation.overlapping_segments.begin(),
				segmentation.overlapping_segments.begin() + handled_count
			);
			libbio_assert(segmentation.overlapping_variants.empty() || !segmentation.overlapping_segments.empty());
			
			// Maintain the sorted order.
			m_variant_writer.merge_output_variants(new_variant_count);
			// Output.
			m_variant_writer.filter_processed_variants_and_output(min_unhandled_ref_pos);
		}
		segmentation.overlap_counter.clean_up_counts(
			segmentation.overlapping_variants.empty()
			? SIZE_MAX
			: segmentation.overlapping_segments.front().alt.position
		);
	}
	
	
	// Entry point.
	void msa_combiner::process_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen)
	{
		m_variant_writer.output_vcf_header();
		m_data_source.process_msa(ref, alt, var_rec_gen, *this); // Calls process_variants().
	}
}
