/*
 * Copyright (c) 2019–2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <range/v3/all.hpp>
#include "msa_data_source.hh"
#include "pairwise_view.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {

	void msa_data_source::push_current_segment()
	{
		// If a segment is mixed, partition it.
		if (segment_type::MIXED == m_current_segment.type || segment_type::MIXED_ALT_STARTS_WITH_GAP == m_current_segment.type)
			split_mixed_segment(m_current_segment, m_segmentation.overlapping_segments);
		else
		{
			check_segment(m_current_segment);
			auto &seg(m_segmentation.overlapping_segments.emplace_back());
			// Move the strings.
			using std::swap;
			swap(m_current_segment, seg);
		}
		
		// seg.alt.position may be SIZE_MAX in case the alt sequence starts with “-”.
		static_assert(std::is_unsigned_v <decltype(m_segmentation.overlapping_segments.front().alt.position)>);
		libbio_assert(
			ranges::is_sorted(
				m_segmentation.overlapping_segments,
				ranges::less(),
				[](auto const &seg){ return 1 + seg.alt.position; }
			)
		);
		
		libbio_assert(m_segmentation_handler);
		m_segmentation_handler->process_variants(m_segmentation);
	}
	
	
	// Handle aligned_character_packs.
	void msa_data_source::handle(aligned_character_pack const &&pack)
	{
		++m_handled_characters;
		
		m_fsm.update_characters(pack);
		parse_msa(pack);
		
		if (m_logs_status && 0 == m_handled_characters % 10000000)
			lb::log_time(std::cerr) << "Handled " << m_handled_characters << " characters…\n";
	}
	
	
	// Handle variant_records.
	void msa_data_source::handle(variant_record &&rec)
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

		m_segmentation.overlap_counter.push_count(rec, m_ploidy);
		m_segmentation.overlapping_variants.emplace_back(std::move(rec));
		
		if (m_logs_status && 0 == m_handled_variants % 100000)
			lb::log_time(std::cerr) << " Handled " << m_handled_variants << " variants…\n";
	}
	
	
	auto msa_data_source::check_gaps_at_start(vector_type const &ref, vector_type const &alt) const -> gap_start_position
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
	void msa_data_source::process_msa(vector_type const &ref, vector_type const &alt, vcf_record_generator &var_rec_gen, segmentation_handler &handler)
	{
		libbio_assert_eq(ref.size(), alt.size());
		
		if (ref.empty()) // implies alt.empty().
			return;
		
		m_segmentation_handler = &handler;
		
		// Iterate over the characters from ref and alt simultaneously and transform them into aligned_strings.
		// At the same time, convert the starting co-ordinate of each variant and store it into a set in the
		// order of the end position relative to alt. When the end co-ordinate is reached in ref and alt, merge the
		// aligned_strings, merge items from the variant set and call the graph handler.
		
		prepare_msa_parser();
		
		libbio_assert(m_segmentation.overlap_counter.empty());
		m_segmentation.overlap_counter.prepare();
		
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
					rsv::single('\0')																	// Mark the end with a special character s.t. the last character will be handled with the sliding window.
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
			[](auto const &var){ return proj_return_type(var.variant.zero_based_pos(), 1); }
		);
		push_current_segment();
		
		m_segmentation_handler = nullptr;
	}
}
