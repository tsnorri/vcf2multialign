/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTIL_HH
#define VCF2MULTIALIGN_UTIL_HH

#include <iostream>


namespace vcf2multialign {
	
	inline void fail(char const *message)
	{
		std::cerr << message << std::endl;
		abort();
	}
	
	
	inline void always_assert(bool check) { if (!check) abort(); }
	
	
	inline void always_assert(bool check, char const *message)
	{
		if (!check)
			fail(message);
	}
	
	
	template <typename t_fn>
	inline void always_assert(bool check, t_fn fn)
	{
		if (!check)
		{
			fn();
			abort();
		}
	}
	
	
	template <typename t_value>
	class copyable_atomic : public std::atomic <t_value>
	{
	public:
		using std::atomic <t_value>::atomic;
		
		copyable_atomic(copyable_atomic const &other)
		{
			this->store(other.load(std::memory_order_acquire), std::memory_order_release);
		}
		
		copyable_atomic &operator=(copyable_atomic const &other) &
		{
			this->store(other.load(std::memory_order_acquire), std::memory_order_release);
			return *this;
		}
	};
	
	
	// From https://stackoverflow.com/a/18940595/856976
	template <typename t_value>
	struct pointer_cmp
	{
		typedef std::true_type is_transparent;
		
		// Use the helper class to turn smart pointers into raw pointers and compare them.
		class helper
		{
		protected:
			t_value const *m_ptr{nullptr};
			
		public:
			helper() = default;
			helper(t_value const &ptr): m_ptr(&ptr) {}
			helper(t_value const *ptr): m_ptr(ptr) {}
			
			template <class ... t_args>
			helper(std::unique_ptr <t_args ...> const &ptr): m_ptr(ptr.get()) {}
			
			template <class ... t_args>
			helper(std::shared_ptr <t_args ...> const &ptr): m_ptr(ptr.get()) {}
			
			bool operator<(helper other) const
			{
				return std::less <t_value const *>()(m_ptr, other.m_ptr);
			}
		};
		
		bool operator()(helper const &&lhs, helper const &&rhs) const
		{
			return lhs < rhs;
		}
	};
	
	
	// Calculate the printed length of a UTF-8 string by checking the first two bits of each byte.
	std::size_t strlen_utf8(std::string const &str);
}

#endif
