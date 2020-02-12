/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_FORWARDER_HH
#define VCF2MULTIALIGN_UTILITY_FORWARDER_HH

#include <cstddef>
#include <stdexcept>
#include <utility>


namespace vcf2multialign {
	
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
