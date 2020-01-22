/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_OVERLAP_COUNTER_HH
#define VCF2MULTIALIGN_COMBINE_MSA_OVERLAP_COUNTER_HH

#include <algorithm>
#include <ostream>
#include <vector>
#include "types.hh"
#include "utility.hh"


namespace vcf2multialign {
	
	class overlap_counter
	{
	public:
		struct overlap_count
		{
			std::size_t		position{};
			std::int32_t	running_sum{};
			std::int8_t		count{};
			
			overlap_count() = default;
			overlap_count(std::size_t const position_, std::int8_t const count_):
				position(position_),
				count(count_)
			{
			}
		};
		
		friend std::ostream &operator<<(std::ostream &, overlap_count const &);
		typedef std::vector <overlap_count>				overlap_count_vector;
		typedef overlap_count_vector::const_iterator	const_iterator;
		typedef const_iterator							iterator;
		
	protected:
		overlap_count_vector				m_overlap_counts;
	
	public:
		bool empty() const { return m_overlap_counts.empty(); }
		void prepare() { m_overlap_counts.emplace_back(0, 0); }
		const_iterator begin() const { return m_overlap_counts.begin(); }
		const_iterator end() const { return m_overlap_counts.end(); }
		overlap_count_vector const &overlap_counts() const { return m_overlap_counts; }
		
		void push_count(variant_record const &var, std::size_t const ploidy);
		void update_running_sums();
		void clean_up_counts(std::size_t const first_unhandled_pos);
		
		inline const_iterator find_initial(std::size_t const seg_alt_pos) const;
		inline const_iterator find_end(const_iterator const overlap_it, std::size_t const seg_alt_end_pos) const;
		
		inline void update_overlap_iterator_if_needed(
			const_iterator &overlap_it,
			std::size_t const seg_alt_pos
		) const;
		
		inline std::size_t initial_count(const_iterator overlap_it) const;
		
		auto max_overlaps_in_range(
			const_iterator overlap_it,
			const_iterator const overlap_end,
			std::size_t const var_pos,
			std::size_t const var_end_pos
		) const -> std::pair <std::int32_t, const_iterator>;
	};
	
	
	inline std::ostream &operator<<(std::ostream &os, overlap_counter::overlap_count const &oc)
	{
		os << "pos: " << oc.position << " rs: " << +oc.running_sum << " count: " << +oc.count;
		return os;
	}
	
	
	void overlap_counter::update_overlap_iterator_if_needed(
		const_iterator &overlap_it,
		std::size_t const seg_alt_pos
	) const
	{
		// Make the overlap range not point to the first position of the segment.
		while (m_overlap_counts.end() != overlap_it && seg_alt_pos == overlap_it->position)
			++overlap_it;
	}
	
	
	auto overlap_counter::find_initial(std::size_t const seg_alt_pos) const -> const_iterator
	{
		auto const overlap_it(
			std::partition_point(
				m_overlap_counts.begin(),
				m_overlap_counts.end(),
				[seg_alt_pos](auto const &oc){ return oc.position < seg_alt_pos; }
			)
		);
		return overlap_it;
	}
	
	
	auto overlap_counter::find_end(const_iterator const overlap_it, std::size_t const seg_alt_end_pos) const -> const_iterator
	{
		auto const overlap_end(
			std::partition_point(
				overlap_it,
				m_overlap_counts.end(),
				[seg_alt_end_pos](auto const &oc){ return oc.position < seg_alt_end_pos; }
			)
		);
		return overlap_end;
	}
	
	
	std::size_t overlap_counter::initial_count(const_iterator overlap_it) const
	{
		auto const initial_overlap_count(m_overlap_counts.begin() == overlap_it ? 0 : (overlap_it - 1)->running_sum);
		return initial_overlap_count;
	}
}

#endif
