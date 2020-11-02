/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include "mnv_combiner.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	void mnv_combiner::handle_variant_description(variant_description &&desc)
	{
		if (!m_previous_desc)
		{
			m_previous_desc.emplace(std::move(desc));
			return;
		}
		
		// Handle the following situation, which is difficult to detect in the MSA parser.
		//
		// POS    REF  ALT
		// p      c    cX
		// p + 1  X    <DEL>
		//
		// Swapping REF and ALT in the first variant and omitting the second should suffice:
		//
		// p      cX   c
		//
		// Check if the current description’s REF matches the ALT in the previous one,
		// as well as the suffix of the previous REF.
		if (
			desc.alt.empty() &&
			desc.position == 1 + m_previous_desc->position &&
			1 + desc.ref.size() == m_previous_desc->alt.size() &&
			m_previous_desc->alt.ends_with(desc.ref)
		)
		{
			using std::swap;
			auto &prev_desc(*m_previous_desc);
			swap(prev_desc.ref, prev_desc.alt);
			
			// Handle and reset.
			m_next_handler->handle_variant_description(std::move(prev_desc));
			m_previous_desc.reset();
		}
		else
		{
			// Store the current description.
			using std::swap;
			swap(*m_previous_desc, desc);
			
			// Handle the previous description.
			m_next_handler->handle_variant_description(std::move(desc));
		}
	}
	
	
	void mnv_combiner::finish()
	{
		if (m_previous_desc)
			m_next_handler->handle_variant_description(std::move(*m_previous_desc));
	}
}
