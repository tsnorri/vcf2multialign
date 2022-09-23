/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include "mnv_combiner.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	
	bool should_rewrite_description(v2m::variant_description const &prev_desc, v2m::variant_description const &desc)
	{
		// Check for
		// POS		REF	ALT
		// p		c	cX
		// p + 1	XY	<DEL>
		
		// Checked by the caller.
		libbio_assert_eq(1 + prev_desc.position, desc.position);
		
		// <DEL>
		if (!desc.alt.empty())
			return false;
		
		// c’s length
		if (1 != prev_desc.ref.size())
			return false;
		
		// cX has c as prefix
		if (!prev_desc.alt.starts_with(prev_desc.ref))
			return false;
		
		// X is a prefix of XY
		std::string_view const alt(prev_desc.alt);
		if (!desc.ref.starts_with(alt.substr(1)))
			return false;
		
		return true;
	}
}


namespace vcf2multialign {
	
	void mnv_combiner::handle_variant_description(variant_description &&desc)
	{
		// Handle the following situation, which is difficult to detect in the MSA parser.
		//
		// POS    REF  ALT
		// p      c    cX
		// p + 1  XY   <DEL>
		//
		// Algorithm
		// =========
		// 1. Sort variants s.t. for variant origin, MSA < VC. (Check the filter class.)
		// 2. Given a position, store the VC variants.
		// 3. When the position changes,
		//	– If the next variant is from VC, pass the stored variants and store the current variant instead.
		//	– If the next variant is from MSA, compare with each of the stored variants.
		//		– If there is no match, pass the stored variants.
		//		– If there is a match, rewrite the stored (VC) variant and continue matching.
		//
		// Rewriting
		// ---------
		// – Check if the suffix of the alternative allele of the previous variant is a prefix of the
		//   reference allele of the current variant (i.e. X is a prefix of XY).
		//		– If not, do not modify the previous variant.
		//		– If yes, swap ALT and REF, flip the genotype, and mark the MSA variant to be discarded.
		
		if (variant_origin::VC == desc.origin)
		{
			// Check if the previous descriptions can be passed to the next handler.
			if (!m_previous_descs.empty() && m_previous_descs.front().position != desc.position)
			{
				for (auto &prev_desc : m_previous_descs)
					m_next_handler->handle_variant_description(std::move(prev_desc));
				
				m_previous_descs.clear();
			}
			
			m_previous_descs.emplace_back(std::move(desc));
		}
		else
		{
			// variant_origin::MSA == desc.origin.
			if (m_previous_descs.empty())
			{
				m_next_handler->handle_variant_description(std::move(desc));
				return;
			}
			
			// Check if we need to try to rewrite the previous variants.
			if (m_previous_descs.front().position + 1 != desc.position)
			{
				for (auto &prev_desc : m_previous_descs)
					m_next_handler->handle_variant_description(std::move(prev_desc));
				
				m_previous_descs.clear();
				m_next_handler->handle_variant_description(std::move(desc));
				return;
			}
			
			// The positions match.
			bool should_skip_current(false);
			for (auto &prev_desc : m_previous_descs)
			{
				if (should_rewrite_description(prev_desc, desc))
				{
					// If the variant caller has found anything similar to the MSA-derived variant,
					// ignore the latter.
					should_skip_current = true;
					prev_desc.genotype.flip();
					
					using std::swap;
					swap(prev_desc.ref, prev_desc.alt);
					++m_combined_variants;

					if (prev_desc.had_alt_eq_to_ref)
					{
						std::cerr << "WARNING: Combining with a variant where ALT matched REF.\n";
						std::cerr << "Prev: " << prev_desc << '\n';
						std::cerr << "Curr: " << desc << '\n';
					}
				}
				
				m_next_handler->handle_variant_description(std::move(prev_desc));
			}
			
			m_previous_descs.clear();
			if (!should_skip_current)
				m_next_handler->handle_variant_description(std::move(desc));
		}
	}
	
	
	void mnv_combiner::finish()
	{
		for (auto &desc : m_previous_descs)
			m_next_handler->handle_variant_description(std::move(desc));
		
		m_previous_descs.clear();
	}
}
