/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include "mnv_combiner.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	bool mnv_combiner::should_combine_with_previous(variant_description const &desc) const
	{
		// Check for
		// POS		REF	ALT
		// p		c	cX
		// p + 1	X	<DEL>
		
		libbio_assert(m_previous_desc);
		auto const &prev_desc(*m_previous_desc);
		
		// <DEL>
		if (!desc.alt.empty())
			return false;
		
		// p and p + 1
		if (desc.position != 1 + prev_desc.position)
			return false;
		
		// c’s length
		if (1 != prev_desc.ref.size())
			return false;
		
		// cX has c as prefix
		if (!prev_desc.alt.starts_with(prev_desc.ref))
			return false;
		
		// cX has X as suffix
		if (!prev_desc.alt.ends_with(desc.ref))
			return false;
		
		return true;
	}
	
	
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
		if (should_combine_with_previous(desc))
		{
			if (1 == m_ploidy)
			{
				auto &prev_desc(*m_previous_desc);
				auto const genotypes{prev_desc.genotype[0] | (std::uint16_t(desc.genotype[0]) << 1)};
				switch (genotypes)
				{
					case 0x0:
					case 0x3:
						// The variants cancel out each other.
						m_previous_desc.reset();
						break;
						
					case 0x1:
						// Both Xs are effective.
						prev_desc.alt += desc.alt;
						break;
						
					case 0x2:
						// cX -> c
						prev_desc = std::move(desc);
						break;
						
					default:
						libbio_fail("Unexpected genotype value");
				}
			}
			else
			{
				// Discard the variants as this is too difficult to handle
				// (unless maybe if the variants were phased).
				*m_previous_desc = std::move(desc);
			}
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
