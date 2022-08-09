/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_RANGE_HH
#define VCF2MULTIALIGN_RANGE_HH

#include <cstddef>


namespace vcf2multialign {
	
	class range
	{
	public:
		typedef std::size_t	size_type;
		
	protected:
		size_type	m_position{};
		size_type	m_length{};
		
	public:
		constexpr range() = default;
		constexpr range(size_type const position, size_t const length):
			m_position(position),
			m_length(length)
		{
		}
		
		constexpr size_type position() const { return m_position; }
		constexpr size_type length() const { return m_length; }
		
		template <typename t_item>
		constexpr void check(std::vector <t_item> const &vec) const;
	};
	
	
	class drop_last
	{
	public:
		typedef std::size_t	size_type;
		
	protected:
		size_type m_count{};
		
	public:
		constexpr drop_last() = default;
		constexpr drop_last(std::size_t const count):
			m_count(count)
		{
		}
		
		constexpr size_type count() const { return m_count; }
	};
	
	
	constexpr inline range
	operator|(range const rng, drop_last const dl)
	{
		auto const len(rng.length());
		return range(rng.position(), len - std::min(len, dl.count()));
	}
	
	
	template <typename t_item>
	constexpr void range::check(std::vector <t_item> const &vec) const
	{
		if (m_position <= m_position + m_length && m_position < vec.size() && m_position + m_length <= vec.size())
			return;
		
		throw std::out_of_range("Given vector indices out of range");
	}
}

#endif
