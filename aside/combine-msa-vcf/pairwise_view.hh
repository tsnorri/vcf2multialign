/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VCF_PAIRWISE_VIEW_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VCF_PAIRWISE_VIEW_HH

#include <range/v3/all.hpp>


namespace vcf2multialign { namespace range {
	
	// I wasn’t able to use ranges::view::sliding with ranges::view::cache1.
	template <typename t_range>
	class pairwise_view : public ranges::view_facade <pairwise_view <t_range>>
	{
		friend ranges::range_access;
		
	public:
		typedef std::true_type													single_pass;
		
		// FIXME: I don’t know why I need to use the latter method to get the value type of the underlying range.
		//typedef	remove_cvref_t <ranges::value_type <t_range>>				underlying_value_type;
		typedef remove_cvref_t <decltype(*std::declval <t_range>().begin())>	underlying_value_type;
		
		typedef	std::pair <underlying_value_type, underlying_value_type>		value_type;
		
	public:
		class cursor;
		
		class sentinel
		{
			friend cursor;
			
		protected:
			ranges::sentinel_t <t_range>	m_sentinel{};
			
		public:
			sentinel() = default;
			
			constexpr explicit sentinel(ranges::sentinel_t <t_range> const sentinel):
				m_sentinel(sentinel)
			{
			}
		};
		
		class cursor
		{
		protected:
			pairwise_view									*m_view{};
			ranges::iterator_t <t_range>					m_it{};
			
		public:
			typedef pairwise_view::value_type				value_type;
			typedef std::true_type							single_pass;
			typedef ranges::range_difference_t <t_range>	difference_type;
			
		public:
			cursor() = default;
			
			constexpr cursor(pairwise_view &view, ranges::iterator_t <t_range> &&it):
				m_view(&view),
				m_it(std::move(it))
			{
			}
			
			inline void next();
			inline value_type const &read() const;
			bool equal(cursor const &other) const { return m_it == other.m_it; }
			bool equal(sentinel const &other) const { return m_it == other.m_sentinel; }
		};
		
	protected:
		t_range		m_range{};
		value_type	m_current_value{};
		bool		m_should_update{};
		
	public:
		pairwise_view() = default;
		
		pairwise_view(t_range &&range):
			m_range(std::forward <t_range>(range))
		{
		}
		
		inline cursor begin_cursor();
		auto end_cursor() { return end_cursor_2(std::bool_constant <ranges::common_range <t_range>>()); }
		
	protected:
		inline void update(ranges::iterator_t <t_range> const &it);
		inline cursor end_cursor_2(std::true_type);
		inline sentinel end_cursor_2(std::false_type);
	};
	
	
	template <typename t_range>
	void pairwise_view <t_range>::cursor::next()
	{
		++m_it;
		m_view->m_should_update = true;
	}
	
	
	template <typename t_range>
	auto pairwise_view <t_range>::cursor::read() const -> value_type const &
	{
		if (m_view->m_should_update)
			m_view->update(m_it);
		return m_view->m_current_value;
	}
	
	
	template <typename t_range>
	void pairwise_view <t_range>::update(ranges::iterator_t <t_range> const &it)
	{
		m_current_value.first = std::move(m_current_value.second);
		m_current_value.second = *it;
		m_should_update = false;
	}
	
	
	template <typename t_range>
	auto pairwise_view <t_range>::begin_cursor() -> cursor
	{
		auto it(ranges::begin(m_range));
		m_current_value.first = *it;
		m_current_value.second = *(++it);
		return cursor(*this, std::move(it));
	}
	
	
	template <typename t_range>
	auto pairwise_view <t_range>::end_cursor_2(std::true_type) -> cursor
	{
		return cursor(ranges::end(m_range));
	}
	
	
	template <typename t_range>
	auto pairwise_view <t_range>::end_cursor_2(std::false_type) -> sentinel
	{
		return sentinel(ranges::end(m_range));
	}
	
	
	struct pairwise_fn
	{
		template <typename t_range>
		auto operator()(t_range &&range) const
		{
			return pairwise_view(std::forward <t_range>(range));
		}
		
		template <typename t_range>
		friend auto operator|(t_range &&range, pairwise_fn const &)
		{
			return pairwise_view(std::forward <t_range>(range));
		}
	};
}}


namespace vcf2multialign { namespace range { namespace view {
	pairwise_fn constexpr pairwise;
}}}

#endif
