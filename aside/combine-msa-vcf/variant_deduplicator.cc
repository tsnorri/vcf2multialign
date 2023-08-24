/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <numeric>
#include <range/v3/all.hpp>
#include "output_handler.hh"
#include "variant_deduplicator.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace {
	class variant_desc_cmp
	{
		static_assert(v2m::variant_origin::MSA < v2m::variant_origin::VC);
		
	public:
		typedef std::tuple <
			std::size_t,
			std::string const &,
			std::string const &,
			v2m::variant_origin
		> tuple_type;
		
	protected:
		static tuple_type to_tuple(v2m::variant_description const &var)
		{
			return {var.position, var.ref, var.alt, var.origin};
		}
		
	public:
		bool operator()(v2m::variant_description const &lhs, v2m::variant_description const &rhs) const
		{
			return to_tuple(lhs) < to_tuple(rhs);
		}
	};
	
	
	typedef std::tuple <
		std::size_t,
		std::string const &,
		std::string const &
	> variant_eq_tuple;
	
	variant_eq_tuple to_tuple_for_eq(v2m::variant_description const &var)
	{
		return {var.position, var.ref, var.alt};
	}
}


namespace vcf2multialign {
	
	void variant_deduplicator::merge_output_variants(std::size_t const partition_point)
	{
		// Merge the partitions of sorted variants.
		variant_desc_cmp cmp;
		std::sort(
			m_output_variants.begin() + partition_point,
			m_output_variants.end(),
			cmp
		);
		std::inplace_merge(
			m_output_variants.begin(),
			m_output_variants.begin() + partition_point,
			m_output_variants.end(),
			cmp
		);
		libbio_assert(
			std::is_sorted(
				m_output_variants.begin(),
				m_output_variants.end(),
				cmp
			)
		);
	}
	
	
	void variant_deduplicator::filter_processed_variants_and_output(std::size_t const min_unhandled_ref_pos)
	{
		// Omit filtering for now, except for checking REF against ALT.
		libbio_assert(std::is_sorted(m_output_variants.begin(), m_output_variants.end(), variant_desc_cmp{}));
		auto var_it(m_output_variants.begin());
		auto const var_end(m_output_variants.end());
		variant_desc_cmp cmp;
		
		if (var_it != var_end)
		{
			// Take the first variant.
			auto *desc(&*var_it);
			if (min_unhandled_ref_pos <= desc->position)
				return;
			++var_it;
			
			// Compare to the next one and pass to the next handler if needed.
			while (var_it != var_end)
			{
				auto &next_desc(*var_it);
				if (min_unhandled_ref_pos <= next_desc.position)
					break;
				
				if (to_tuple_for_eq(*desc) != to_tuple_for_eq(next_desc))
					m_next_handler->handle_variant_description(std::move(*desc));
				
				// MSA < VC, so we can pass whichever variant is stored last.
				desc = &next_desc;
				++var_it;
			}
			
			m_next_handler->handle_variant_description(std::move(*desc));
			m_output_variants.erase(m_output_variants.begin(), var_it);
		}
	}
	
	
	void variant_deduplicator::finish()
	{
		filter_processed_variants_and_output(SIZE_MAX);
		m_next_handler->finish();
	}
}
