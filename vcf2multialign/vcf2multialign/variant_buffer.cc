/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <vcf2multialign/variant_buffer.hh>


namespace vcf2multialign {

	auto variant_buffer::variant_range() -> std::pair <variant_list::iterator, variant_list::iterator>
	{
		auto it(m_variant_list.begin());
		return std::make_pair(it, it + m_list_ptr);
	}


	void variant_buffer::fill_buffer()
	{
		auto const size(m_variant_list.size());

		// If there is an item that has not been handled before, move it to the beginning of the list.
		if (m_list_ptr)
		{
			using std::swap;
			swap(m_variant_list[m_list_ptr], m_variant_list[0]);
			m_list_ptr = 1;
		}

		while (true)
		{
			// Copy a model variant if needed.
			if (! (m_list_ptr < size))
				m_variant_list.push_back(m_model_variant);

			auto &current_variant(m_variant_list[m_list_ptr]);

			// Read the line.
			if (!m_reader->get_line(current_variant.line))
				break;

			// Parse the variant.
			m_reader->get_next_variant(current_variant.var, current_variant.line);

			// Check if the new variant has the same position as the previous ones.
			if (m_variant_list[0].var.pos() != current_variant.var.pos()) 
				break;

			++m_list_ptr;
		}

		// Sort so that the longest REF is in the beginning of the list.
		auto const begin(m_variant_list.begin());
		std::sort(begin, begin + m_list_ptr, [](buffered_variant &lhs, buffered_variant &rhs) {
			return lhs.var.ref().size() > rhs.var.ref().size();
		});
	}
}
