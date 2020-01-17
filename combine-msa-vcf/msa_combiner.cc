/*
 * Copyright (c) 2019–2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/variant_printer.hh>
#include <range/v3/all.hpp>
#include <vcf2multialign/variant_format.hh>
#include "algorithms.hh"
#include "msa_combiner.hh"
#include "pairwise_view.hh"
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
}


namespace vcf2multialign {

	void msa_combiner::push_current_segment()
	{
		// If a segment is mixed, partition it.
		if (segment_type::MIXED == m_current_segment.type)
			split_mixed_segment(m_current_segment, m_overlapping_segments);
		else
		{
			check_segment(m_current_segment);
			auto &seg(m_overlapping_segments.emplace_back());
			// Move the strings.
			using std::swap;
			swap(m_current_segment, seg);
		}
		
		// seg.alt.position may be SIZE_MAX in case the alt sequence starts with “-”.
		static_assert(std::is_unsigned_v <decltype(m_overlapping_segments.front().alt.position)>);
		libbio_assert(
			ranges::is_sorted(
				m_overlapping_segments,
				ranges::less(),
				[](auto const &seg){ return 1 + seg.alt.position; }
			)
		);
		
		process_variants();
	}
	
	
	// Handle aligned_character_packs.
	void msa_combiner::handle(aligned_character_pack const &&pack)
	{
		++m_handled_characters;
		
		m_fsm.update_characters(pack);
		parse_msa(pack);
		
		if (m_logs_status && 0 == m_handled_characters % 1000)
			lb::log_time(std::cerr) << " Handled " << m_handled_characters << " characters…\n";
	}
	
	
	// Handle variant_records.
	void msa_combiner::handle(variant_record &&rec)
	{
		++m_handled_variants;
		
		// Transform the co-ordinate value in rec.
		auto const rec_alt_pos(rec.aligned_position - 1); // The position relative to alt; rec.aligned_position is not aligned yet.
		auto const rec_len(rec.variant.ref().size());
		// The current or previous segment needs to be updated before this can be done,
		// so handle() for variant_record cannot be called first.
		libbio_always_assert(m_current_segment.alt.position <= rec_alt_pos);
		libbio_assert(segment_type::MATCH == m_current_segment.type || !m_current_segment.alt.string.empty());	// The variant may not be aligned to a deleted position.
		libbio_assert_lte(m_current_segment.alt.position, rec_alt_pos);
		auto const rec_pos_in_seg(rec_alt_pos - m_current_segment.alt.position);
		rec.aligned_position = m_current_segment.aligned_position + rec_pos_in_seg;

		m_overlap_counter.push_count(rec, m_ploidy);
		m_overlapping_variants.emplace_back(std::move(rec));
		
		if (m_logs_status && 0 == m_handled_variants % 1000)
			lb::log_time(std::cerr) << " Handled " << m_handled_variants << " variants…\n";
	}
	
	
	void msa_combiner::handle_one_segment_msa(
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
					m_variant_writer.emplace_back(
						variant_description(
							seg.ref.position + begin_idx,
							ref_sv.substr(begin_idx, end_idx - begin_idx),
							alt_sv.substr(begin_idx, end_idx - begin_idx),
							m_ploidy,
							overlap_count,
							variant_origin::MSA
						)
					);
					
					begin_idx = end_idx;
					overlap_count = oc.running_sum;
				}
				
				if (begin_idx < alt_sv.size())
				{
					m_variant_writer.emplace_back(
						variant_description(
							seg.ref.position + begin_idx,
							ref_sv.substr(begin_idx),
							alt_sv.substr(begin_idx),
							m_ploidy,
							overlap_count,
							variant_origin::MSA
						)
					);
				}
				break;
			}
				
			case segment_type::DELETION:
			{
				libbio_assert_eq(overlap_it, overlap_end);
				m_variant_writer.emplace_back(
					variant_description(
						seg.ref.position,
						seg.ref.string,
						"",
						m_ploidy,
						overlap_count,
						variant_origin::MSA
					)
				);
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
				m_variant_writer.emplace_back(
					variant_description(
						seg.ref.position,
						seg.ref.string,
						seg.alt.string,
						m_ploidy,
						max_overlaps,
						variant_origin::MSA
					)
				);
				break;
			}
				
			case segment_type::MIXED:
				throw std::runtime_error("Unexpected segment type");
		}
	}
	
	
	std::size_t msa_combiner::process_variants_msa(
		std::size_t const max_alt_pos,
		aligned_segment_vector::const_iterator seg_it,
		aligned_segment_vector::const_iterator const seg_end
	)
	{
		if (seg_it == seg_end) return 0;
		
		// Determine the overlap count range for the first segment.
		auto const &first_seg(*seg_it);
		auto const first_seg_alt_pos(first_seg.alt.position);
		auto const first_seg_alt_end_pos(first_seg.alt_end());
		auto overlap_it(m_overlap_counter.find_initial(first_seg_alt_pos));
		
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
			
			m_overlap_counter.update_overlap_iterator_if_needed(overlap_it, seg_alt_pos);
			auto const initial_overlap_count(m_overlap_counter.initial_count(overlap_it));
			auto const overlap_end(m_overlap_counter.find_end(overlap_it, seg_alt_end_pos));
			handle_one_segment_msa(seg, initial_overlap_count, overlap_it, overlap_end);
			
			overlap_it = overlap_end;
			++handled_count;
			++seg_it;
		}
		
		return handled_count;
	}
	
	
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
		//				   YY   YY
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
		// Var2	  TTC		  TCC	<- Variants need to be rewritten like this.
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
		auto const &last_seg(*(aligned_segment_end - 1));
		libbio_assert_neq(first_seg.type, segment_type::DELETION);
		libbio_assert_neq(last_seg.type, segment_type::DELETION);
		auto const var_pos_seg_relative(segment_relative_variant_position(first_seg, var));
		
		// Handle the variant.
		{
			auto const range(ranges::subrange(aligned_segment_begin, aligned_segment_end));
			std::size_t alt_idx{};
			for (auto const &var_alt : var.variant.alts())
			{
				++alt_idx;
				auto left_pad(var_pos_seg_relative);
				std::size_t var_ref_characters_remaining(var.size); // Number of remaining variant reference characters, i.e. in seg.alt.
				std::size_t var_alt_characters_remaining(var_alt.alt.size());
				std::size_t total_alt_characters_consumed{};
				auto const [gt_count, ploidy] = (count_set_genotype_values(var.variant, alt_idx));
				libbio_always_assert_eq_msg(ploidy, m_ploidy, "Line ", var.variant.lineno(), ": expected the sample ploidy to match the passed value, got ", ploidy, '.');
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
			}
		}
	}
	
	
	// Handle an overlapping part.
	void msa_combiner::process_variants()
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
		// reported by the VC. Even though this should only be called with a group of overlapping variants,
		// the fact that two VC-provided variants do or do not overlap does not affect anything.
		
		m_overlap_counter.update_running_sums();
		auto const initial_variant_count(m_variant_writer.size());
		if (m_overlapping_variants.empty())
		{
			// Handle MSA only.
			auto const end(m_overlapping_segments.end());
			auto const handled_count(process_variants_msa(SIZE_MAX, m_overlapping_segments.cbegin(), end));
			m_overlapping_segments.erase(m_overlapping_segments.begin(), m_overlapping_segments.begin() + handled_count);
			
			m_variant_writer.merge_output_variants(initial_variant_count);
			m_variant_writer.filter_processed_variants_and_output(SIZE_MAX);
		}
		else
		{
			libbio_assert(!m_overlapping_segments.empty());
			libbio_assert_neq(SIZE_MAX, m_overlapping_segments.front().alt.position);
			
			std::size_t min_unhandled_alt_pos(SIZE_MAX);
			std::size_t min_unhandled_ref_pos(SIZE_MAX);
			
			auto var_it(m_overlapping_variants.begin());
			auto seg_it(m_overlapping_segments.begin());
			auto overlap_it(m_overlap_counter.begin());
			auto const var_end(m_overlapping_variants.end());
			auto const seg_end(m_overlapping_segments.end());
			auto const overlap_end(m_overlap_counter.end());
			
			while (var_it != var_end)
			{
				auto const &var(*var_it);
				auto const var_pos(var.variant.zero_based_pos());
				auto const var_end_pos(var_pos + var.size);
				
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
				auto const res(m_overlap_counter.max_overlaps_in_range(overlap_it, overlap_end, var_pos, var_end_pos));
				auto const max_overlaps(res.first);
				overlap_it = res.second;
				process_variant_in_range(var, seg_it, seg_current_end, max_overlaps);
				++var_it;
			}
			
			m_overlapping_variants.erase(
				std::remove_if(
					m_overlapping_variants.begin(),
					m_overlapping_variants.end(),
					[](auto const &var){
						return !var.is_skipped;
					}
				),
				m_overlapping_variants.end()
			);
			for (auto &var : m_overlapping_variants)
				var.is_skipped = false;
			
			m_variant_writer.merge_output_variants(initial_variant_count);
			
			auto const updated_variant_count(m_variant_writer.size());
			auto const handled_count(process_variants_msa(min_unhandled_alt_pos, m_overlapping_segments.begin(), m_overlapping_segments.end()));
			m_overlapping_segments.erase(m_overlapping_segments.begin(), m_overlapping_segments.begin() + handled_count);
			libbio_assert(m_overlapping_variants.empty() || !m_overlapping_segments.empty());
			
			m_variant_writer.merge_output_variants(updated_variant_count);
			m_variant_writer.filter_processed_variants_and_output(min_unhandled_ref_pos);
		}
		m_overlap_counter.clean_up_counts(
			m_overlapping_variants.empty()
			? SIZE_MAX
			: m_overlapping_segments.front().alt.position
		);
	}
	
	
	auto msa_combiner::check_gaps_at_start(vector_type const &ref, vector_type const &alt) const -> gap_start_position
	{
		auto const ref_begin(ref.begin());
		auto const alt_begin(alt.begin());
		auto const ref_end(ref.end());
		auto const alt_end(alt.end());
		auto const ref_it(std::find_if_not(ref_begin, ref_end, [](auto const c){ return '-' == c; }));
		auto const alt_it(std::find_if_not(alt_begin, alt_end, [](auto const c){ return '-' == c; }));
		std::size_t const ref_dist(ref_it == ref_end ? 0 : std::distance(ref_begin, ref_it));
		std::size_t const alt_dist(alt_it == alt_end ? 0 : std::distance(alt_begin, alt_it));
		
		if (ref_dist <= alt_dist)
			return gap_start_position(ref_dist, 0, '\0', false);
		else // alt_dist < ref_dist
			return gap_start_position(alt_dist, ref_dist - alt_dist, *ref_it, true);
	}
	
	
	// Entry point.
	void msa_combiner::process_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen)
	{
		libbio_assert_eq(ref.size(), alt.size());
		
		if (ref.empty()) // implies alt.empty().
			return;
		
		// Iterate over the characters from ref and alt simultaneously and transform them into aligned_strings.
		// At the same time, convert the starting co-ordinate of each variant and store it into a set in the
		// order of the end position relative to alt. When the end co-ordinate is reached in ref and alt, merge the
		// aligned_strings, merge items from the variant set and call the graph handler.
		
		prepare_msa_parser();
		
		libbio_assert(m_overlap_counter.empty());
		m_overlap_counter.prepare();
		
		// First, check if ref starts with a sequence of gap characters and handle it.
		auto const sp_info(check_gaps_at_start(ref, alt));
		m_current_segment.aligned_position = sp_info.aligned_start;
		
		position_counter ref_cnt, alt_cnt;
		
		auto const drop_count(sp_info.aligned_start + sp_info.ref_pad_diff + sp_info.ref_chr_needs_left_align);
		auto lrange(
			rsv::zip(
				rsv::iota(sp_info.aligned_start),
				rsv::concat( // Ref
					rsv::repeat_n(sp_info.first_ref_chr_after_pad, sp_info.ref_chr_needs_left_align),	// If alt has less gap characters in the beginning than ref, align ref’s first
					rsv::repeat_n('-', sp_info.ref_pad_diff),											// character to the left. (This is consistent with VCF 4.3 section 1.6.1.4.)
					ref | rsv::drop_exactly(drop_count),												// Drop also the first non-gap character if it was in the rsv::repeat_n above.
					rsv::single('-')																	// Add one gap to the end s.t. the last character will be handled with the sliding window.
				)
				| rsv::transform([&ref_cnt](auto const c) { return ref_cnt.add_chr(c); }),	// Wrap inside sequence_character.
				rsv::concat( // Alt / ad-hoc ref
					alt | rsv::drop_exactly(sp_info.aligned_start),
					rsv::single('-')
				)
				| rsv::transform([&alt_cnt](auto const c) { return alt_cnt.add_chr(c); })
			)
			| v2m::range::view::pairwise	// Using the sliding window from ranges caused problems.
			| rsv::transform(
				[](auto const &t) { return aligned_character_pack(				// Simplify access by wrapping into aligned_character_pack.
					std::get <0>(t.first),				// Position
					std::get <1>(t.first),				// sequence_character
					std::get <2>(t.first),				// sequence_character
					std::get <1>(t.second).character,	// next character in ref
					std::get <2>(t.second).character);	// next character in alt
				})
				// No need to move the structs since they only contain primitive types.
		);
		auto rrange(
			rsv::generate([&var_rec_gen](){ return var_rec_gen.next_variant(); })
			| rsv::take_while([](auto const &var){ return 0 != var.aligned_position; })		// Call next_variant() while the returned positions are valid.
			| rsv::move																		// Cast to rvalue references since the variants contain allocated strings.
		);
		
		// Merge and pass to this->handle(). Sort the items s.t. lrsv comes first, except if the alt character is a gap.
		// This causes the variant’s position to always be in m_current_segment.
		m_variant_writer.output_vcf_header();
		if (m_logs_status)
			lb::log_time(std::cerr) << "Creating segments and merging…\n";
		forwarder fwd(*this);
		typedef std::tuple <std::size_t, std::uint8_t> proj_return_type;
		ranges::merge(
			lrange,
			rrange,
			fwd,
			ranges::less(),
			[](auto const &pck){ return proj_return_type(pck.alt.position, ('-' == pck.alt.character ? 2 : 0)); },
			[](auto const &var){ return proj_return_type(var.aligned_position - 1, 1); }	// Not yet aligned.
		);
		push_current_segment();
	}
}
