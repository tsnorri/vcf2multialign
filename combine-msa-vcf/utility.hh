/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_UTILITY_HH
#define VCF2MULTIALIGN_COMBINE_MSA_UTILITY_HH

#include <libbio/assert.hh>
#include <libbio/vcf/variant.hh>
#include <range/v3/all.hpp>
#include "algorithms.hh"
#include "types.hh"


namespace vcf2multialign {
	
	// From C++20.
	template <typename t_type>
	struct remove_cvref
	{
		typedef std::remove_cv_t <std::remove_reference_t <t_type>> type;
	};
	
	template <typename t_type>
	using remove_cvref_t = typename remove_cvref <t_type>::type;
	
	
	std::pair <std::uint16_t, std::uint16_t> count_set_genotype_values(libbio::variant const &var, std::uint16_t const alt_idx);
	void split_mixed_segment(aligned_segment &src, aligned_segment_vector &dst);
	
	
	inline void check_segment(aligned_segment const &seg)
	{
		namespace rsv	= ranges::view;
		
		// Check that for MISMATCH, every (aligned) character does not match.
		libbio_assert(
			segment_type::MISMATCH != seg.type ||
			ranges::all_of(
				rsv::zip(seg.ref.string, seg.alt.string),
				[](auto const &tup){
					return std::get <0>(tup) != std::get <1>(tup);
				}
			)
		);
	}
	
	
	// Find the segments that overlap with the given variant.
	// The range from seg_it to seg_end needs to consist of consequtive segments sorted by their position.
	template <typename t_iterator, typename t_variant>
	inline auto find_overlapping_segment_range(
		t_iterator const seg_it,
		t_iterator const seg_end,
		t_variant const &var
	) -> std::pair <t_iterator, t_iterator>
	{
		struct {
			bool operator()(t_variant const &var, aligned_segment const &seg) const
			{
				// lt. if the variant is located before this segment.
				auto const var_pos(var.variant.zero_based_pos());
				auto const var_end_pos(var_pos + var.size);
				return var_end_pos <= seg.alt.position;
			}
			
			bool operator()(aligned_segment const &seg, t_variant const &var) const
			{
				// lt. if the variant start is located after this segment.
				auto const var_pos(var.variant.zero_based_pos());
				auto const seg_end_pos(seg.alt_end());
				return seg_end_pos <= var_pos;
			}
		} seg_cmp;
		
		return std::equal_range(seg_it, seg_end, var, seg_cmp);
	}
	
	
	// May be used in place of an output iterator to forward arguments of operator=.
	template <typename t_dst>
	class forwarder
	{
	public:
		typedef std::ptrdiff_t difference_type; // Needed for range-v3 only.
		
	protected:
		t_dst	*m_dst{};
	
	public:
		// ranges::semiregular<T> (in range-v3/include/concepts/concepts.hpp) requires copyable and default_constructible.
		// m_dst is required, though, so throw in case the default constructor is somehow called.
		forwarder()
		{
			throw std::runtime_error("forwarder’s default constructor should not be called.");
		}
		
		forwarder(t_dst &dst):
			m_dst(&dst)
		{
		}
	
		forwarder &operator++() { return *this; }				// Return *this.
		forwarder &operator++(int) { return *this; }			// Return *this.
		forwarder &operator*() { return *this; }				// Return *this.
		forwarder const &operator*() const { return *this; }	// Return *this.
	
		template <typename t_arg>
		forwarder &operator=(t_arg &&arg) { m_dst->handle(std::forward <t_arg>(arg)); return *this; }
	};
}

#endif