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
	boost::bimaps::set_of <size_t>,		// lineno
	boost::bimaps::list_of <size_t>		// count
> conflict_count_map;


namespace {
	
	struct var_info
	{
		size_t pos;
		size_t lineno;
		
		var_info(size_t const pos_, size_t const lineno_):
			pos(pos_),
			lineno(lineno_)
		{
		}
	};

	
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

	size_t check_overlapping_non_nested_variants(
		vcf_reader &reader,
		variant_set /* out */ &skipped_variants,
		error_logger &error_logger
	)
	{
		typedef boost::bimap <
			boost::bimaps::multiset_of <size_t>,
			boost::bimaps::multiset_of <size_t>
		> overlap_map;
		
		size_t last_position(0);
		std::multimap <size_t, var_info> end_positions; // end -> pos & lineno
		conflict_count_map conflict_counts;
		overlap_map bad_overlaps;
		size_t i(0);
		size_t conflict_count(0);
		
		reader.reset();
		reader.set_parsed_fields(vcf_field::REF);
		bool should_continue(false);
		do {
			reader.fill_buffer();
			should_continue = reader.parse(
				[&last_position, &end_positions, &conflict_counts, &bad_overlaps, &i, &conflict_count]
				(v2m::transient_variant const &var)
				-> bool
			{
				// Verify that the positions are in increasing order.
				auto const pos(var.zero_based_pos());

				if (! (last_position <= pos))
					throw std::runtime_error("Positions not in increasing order");

				auto const end(pos + var.ref().size());
				auto const var_lineno(var.lineno());

				// Try to find an end position that is greater than var's position.
				auto it(end_positions.upper_bound(pos));
				auto const end_it(end_positions.cend());

				if (end_it == it)
				{
					end_positions.emplace(
						std::piecewise_construct,
						std::forward_as_tuple(end),
						std::forward_as_tuple(pos, var_lineno)
					);
					goto loop_end;
				}
			
				// Check that the found position is not within var's range.
				// If it is, continue checking succeeding positions.
				do
				{
					// Proper nesting.
					if (end <= it->first)
						break;
				
					// Check if the potentially conflicting variant is in fact inside this one.
					auto const other_lineno(it->second.lineno);
					auto const other_pos(it->second.pos);
					if (pos == other_pos)
						goto loop_end_2;
				
					++conflict_count;
					std::cerr << "Variant on line " << var_lineno << " conflicts with line " << other_lineno << "." << std::endl;
				
					{
						auto const res(bad_overlaps.insert(overlap_map::value_type(other_lineno, var_lineno)));
						if (false == res.second)
							throw std::runtime_error("Unable to insert");
					}

					++conflict_counts.left[other_lineno];
					++conflict_counts.left[var_lineno];
				
				loop_end_2:
					++it;
				} while (end_it != it);
			
				// Add the end position.
				end_positions.emplace(
					std::piecewise_construct,
					std::forward_as_tuple(end),
					std::forward_as_tuple(pos, var_lineno)
				);

			loop_end:
				last_position = pos;
			
				++i;
				if (0 == i % 100000)
					std::cerr << "Handled " << i << " variants…" << std::endl;
			
				return true;
			});
		} while (should_continue);
		
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
