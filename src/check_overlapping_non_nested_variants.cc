/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <iostream>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>


namespace vcf2multialign {
	
	size_t check_overlapping_non_nested_variants(vcf_reader &reader, variant_set /* out */ &skipped_variants)
	{
		typedef boost::bimap <
			boost::bimaps::multiset_of <size_t>,
			boost::bimaps::multiset_of <size_t>
		> overlap_map;
		typedef boost::bimap <
			boost::bimaps::set_of <size_t>,
			boost::bimaps::list_of <size_t>
		> conflict_count_map;
		
		std::cerr << "Checking overlapping variants…" << std::endl;
		size_t last_position(0);
		std::map <size_t, size_t> end_positions;
		conflict_count_map conflict_counts;
		overlap_map bad_overlaps;
		size_t i(0);
		size_t conflict_count(0);
		
		reader.reset();
		reader.set_parsed_fields(vcf_field::REF);
		variant var(reader.sample_count());
		while (reader.get_next_variant(var))
		{
			// Verify that the positions are in increasing order.
			auto const pos(var.zero_based_pos());
			if (! (last_position <= pos))
				throw std::runtime_error("Positions not in increasing order");
			
			// Try to find an end position that is greater than var's position.
			auto const end(pos + var.ref().size());
			auto it(end_positions.upper_bound(pos));
			auto const end_it(end_positions.cend());
			auto const var_lineno(var.lineno());

			if (end_it == it)
			{
				end_positions.emplace(end, var_lineno);
				goto loop_end;
			}
			
			// Check that the found position is not within var's range.
			// If it is, continue checking succeeding positions.
			do
			{
				// Proper nesting.
				if (end <= it->first)
					break;
				
				++conflict_count;
				std::cerr << "Variant on line " << var_lineno << " conflicts with line " << it->second << "." << std::endl;
				
				auto const res(bad_overlaps.insert(overlap_map::value_type(it->second, var_lineno)));
				if (false == res.second)
					throw std::runtime_error("Unable to insert");

				++conflict_counts.left[it->second];
				++conflict_counts.left[var_lineno];
				++it;
			} while (end_it != it);
			
			// Add the end position.
			end_positions.emplace(end, var_lineno);

		loop_end:
			++i;
			if (0 == i % 100000)
				std::cerr << "Handled " << i << " variants…" << std::endl;
		}
		
		
		// Remove conflicting variants starting from the one with the highest score.
		// FIXME: currently too many variants are removed.
		conflict_counts.right.sort();
		for (auto const &kv : boost::adaptors::reverse(conflict_counts.right))
		{
			auto const count(kv.first);
			auto const var_lineno(kv.second);
			
			// Check if the candidate variant is still listed.
			{
				auto const left_range(bad_overlaps.left.equal_range(var_lineno));
				if (! (bad_overlaps.left.end() == left_range.first || left_range.first == left_range.second))
				{
					skipped_variants.insert(var_lineno);
					bad_overlaps.left.erase(left_range.first, left_range.second);
				}
			}
			
			{
				auto const right_range(bad_overlaps.right.equal_range(var_lineno));
				if (! (bad_overlaps.right.end() == right_range.first || right_range.first == right_range.second))
				{
					skipped_variants.insert(var_lineno); // May be done b.c. skipped_variants is a set.
					bad_overlaps.right.erase(right_range.first, right_range.second);
				}
			}
			
			if (bad_overlaps.size() == 0)
				break;
		}
		
		if (bad_overlaps.size() != 0)
			throw std::runtime_error("Unable to remove all conflicting variants");
		
		return conflict_count;
	}
}
