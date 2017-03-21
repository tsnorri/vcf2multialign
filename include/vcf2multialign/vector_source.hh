/*
 * Copyright (c) 2016-2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VECTOR_SOURCE_HH
#define VCF2MULTIALIGN_VECTOR_SOURCE_HH

#include <cassert>
#include <iostream>
#include <memory>
#include <mutex>
#include <vector>


namespace vcf2multialign {
	
	template <typename t_vector>
	class vector_source
	{
	public:
		typedef t_vector	vector_type;
		
	protected:
		std::vector <std::unique_ptr <vector_type>>	m_store;
		std::mutex									m_mutex;
		std::size_t									m_in_use{0};
		bool										m_allow_resize{true};
		
	protected:
		void resize(std::size_t const size);
		
	public:
		vector_source(std::size_t size = 0, bool allow_resize = true):
			m_store(size)
		{
			resize(size);
			m_allow_resize = allow_resize;
		}
		
		vector_source(vector_source &&other):
			m_store(std::move(other.m_store)),
			m_allow_resize(other.m_allow_resize)
		{
		}
		
		bool operator==(vector_source const &other) const
		{
			return this == &other;
		}
		
		vector_source &operator=(vector_source &&other) &
		{
			if (*this != other)
			{
				m_store = std::move(other.m_store);
				m_allow_resize = other.m_allow_resize;
			}
			return *this;
		}
		
		void get_vector(std::unique_ptr <vector_type> &target_ptr);
		void put_vector(std::unique_ptr <vector_type> &source_ptr);
	};


	template <typename t_vector>
	void vector_source <t_vector>::resize(std::size_t const size)
	{
		if (!m_allow_resize)
			throw std::runtime_error("Trying to allocate more vectors than allowed");
		
		// Fill from the beginning so that m_in_use slots in the end remain empty.
		m_store.resize(size);
		for (decltype(m_in_use) i(0); i < size - m_in_use; ++i)
			m_store[i].reset(new vector_type);
	}


	template <typename t_vector>
	void vector_source <t_vector>::get_vector(std::unique_ptr <vector_type> &target_ptr)
	{
		assert(nullptr == target_ptr.get());
		
		std::lock_guard <std::mutex> lock_guard(m_mutex);
		auto total(m_store.size());
		
		// Check if there are any vectors left.
		if (m_in_use == total)
		{
			auto new_size(2 * total);
			if (!new_size)
				new_size = 1;
				
			resize(new_size);
			if (m_in_use)
				assert(m_store[m_in_use - 1].get());
			total = new_size;
		}
		
		auto &ptr(m_store[total - m_in_use - 1]);
		target_ptr.swap(ptr);
		++m_in_use;
	}


	template <typename t_vector>
	void vector_source <t_vector>::put_vector(std::unique_ptr <vector_type> &source_ptr)
	{
		assert(source_ptr.get());
		
		std::lock_guard <std::mutex> lock_guard(m_mutex);
		auto const total(m_store.size());
		assert(total);
		assert(m_in_use);
		
		auto &ptr(m_store[total - m_in_use]);
		assert(nullptr == ptr.get());
		source_ptr.swap(ptr);
		--m_in_use;
	}
}

#endif
