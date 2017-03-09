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


namespace v2m = vcf2multialign;


typedef boost::bimap <
	boost::bimaps::multiset_of <size_t>,
	boost::bimaps::multiset_of <size_t>
> overlap_map;


typedef boost::bimap <
	boost::bimaps::set_of <size_t>,
	boost::bimaps::list_of <size_t>
> conflict_count_map;


namespace {
	
	template <typename t_map>
	void check_overlap(
		t_map &bad_overlap_side,
		conflict_count_map &conflict_counts,
		v2m::variant_set &skipped_variants,
		size_t const var_lineno,
		v2m::error_logger &error_logger
	)
	{
		auto const range(bad_overlap_side.equal_range(var_lineno));
		if (! (bad_overlap_side.end() == range.first || range.first == range.second))
		{
			if (error_logger.is_logging_errors())
			{
				for (auto it(range.first); it != range.second; ++it)
					error_logger.log_conflicting_variants(var_lineno, it->second);
			}
			
			// Update conflict counts.
			for (auto it(range.first); it != range.second; ++it)
			{
				auto c_it(conflict_counts.left.find(it->second));
				if (conflict_counts.left.end() == c_it)
					throw std::runtime_error("Unable to find conflict count for variant");
				
				auto &val(c_it->second);
				--val;
				
				if (0 == val)
					conflict_counts.left.erase(c_it);
				
				// In case 0 == val, bad_overlaps need not be updated b.c. the entries have
				// already been erased as part of handling previous overlapping variants.
			}
			
			skipped_variants.insert(var_lineno); // May be done without checking b.c. skipped_variants is a set.
			bad_overlap_side.erase(range.first, range.second);
		}
	}
}


namespace vcf2multialign {
	
	size_t check_overlapping_non_nested_variants(vcf_reader &reader, variant_set /* out */ &skipped_variants, error_logger &error_logger)
	{
		typedef boost::bimap <
			boost::bimaps::multiset_of <size_t>,
			boost::bimaps::multiset_of <size_t>
		> overlap_map;
		
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
		while (!conflict_counts.empty())
		{
			conflict_counts.right.sort();

			auto const it(conflict_counts.right.rbegin());
			auto const count(it->first);
			auto const var_lineno(it->second);
			
			// Check if the candidate variant is still listed.
			check_overlap(bad_overlaps.left, conflict_counts, skipped_variants, var_lineno, error_logger);
			check_overlap(bad_overlaps.right, conflict_counts, skipped_variants, var_lineno, error_logger);
			conflict_counts.left.erase(var_lineno);
		}
		
		if (bad_overlaps.size() != 0)
			throw std::runtime_error("Unable to remove all conflicting variants");
		
		return conflict_count;
	}
}
