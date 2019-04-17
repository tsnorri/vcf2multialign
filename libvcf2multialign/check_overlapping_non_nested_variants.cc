/*
 * Copyright (c) 2017-2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <iostream>
#include <libbio/assert.hh>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


typedef boost::bimap <
	boost::bimaps::multiset_of <std::size_t>,
	boost::bimaps::multiset_of <std::size_t>
> overlap_map;


typedef boost::bimap <
	boost::bimaps::set_of <std::size_t>,		// lineno
	boost::bimaps::list_of <std::size_t>		// count
> conflict_count_map;


namespace {
	
	struct var_info
	{
		std::size_t pos;
		std::size_t lineno;
		
		var_info(std::size_t const pos_, std::size_t const lineno_):
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
		std::size_t const var_lineno,
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
				libbio_always_assert_msg(conflict_counts.left.end() != c_it, "Unable to find conflict count for variant");
				
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

	bool can_handle_variant_alts(lb::transient_variant const &var)
	{
		// Keep.
		for (auto const &alt : var.alts())
		{
			auto const svt(alt.alt_sv_type);
			switch (svt)
			{
				// These structural variant types are currently handled.
				case lb::sv_type::NONE:
				case lb::sv_type::DEL:
				case lb::sv_type::DEL_ME:
					return true;
					
				case lb::sv_type::INS:
				case lb::sv_type::DUP:
				case lb::sv_type::INV:
				case lb::sv_type::CNV:
				case lb::sv_type::DUP_TANDEM:
				case lb::sv_type::INS_ME:
				case lb::sv_type::UNKNOWN:
					break;
			}
		}
		
		return false;
	}
	
	
	bool can_handle_variant_alt(lb::variant_alt_base const &alt)
	{
		auto const svt(alt.alt_sv_type);
		switch (svt)
		{
			case lb::sv_type::NONE:
			case lb::sv_type::DEL:
			case lb::sv_type::DEL_ME:
				return true;
				
			case lb::sv_type::INS:
			case lb::sv_type::DUP:
			case lb::sv_type::INV:
			case lb::sv_type::CNV:
			case lb::sv_type::DUP_TANDEM:
			case lb::sv_type::INS_ME:
			case lb::sv_type::UNKNOWN:
				return false;
		}
	}
	

	std::size_t check_overlapping_non_nested_variants(
		lb::vcf_reader &reader,
		std::string const &chromosome_name,
		variant_set /* out */ &skipped_variants,
		error_logger &error_logger
	)
	{
		typedef boost::bimap <
			boost::bimaps::multiset_of <std::size_t>,
			boost::bimaps::multiset_of <std::size_t>
		> overlap_map;
		
		std::size_t last_position(0);
		std::multimap <std::size_t, var_info> end_positions; // end -> pos & lineno
		conflict_count_map conflict_counts;
		overlap_map bad_overlaps;
		std::size_t i(0);
		std::size_t conflict_count(0);
		
		reader.reset();
		reader.set_parsed_fields(lb::vcf_field::ALT);
		bool should_continue(false);
		do {
			reader.fill_buffer();
			should_continue = reader.parse(
				[
					&chromosome_name,
					&skipped_variants,
					&error_logger,
					&last_position,
					&end_positions,
					&conflict_counts,
					&bad_overlaps,
					&i,
					&conflict_count
				]
				(lb::transient_variant const &var)
				-> bool
			{
				// Verify that the positions are in increasing order.
				auto const pos(var.zero_based_pos());
				libbio_always_assert_msg(last_position <= pos, "Positions not in increasing order");
				
				if (chromosome_name.size() && var.chrom_id() != chromosome_name)
					goto loop_end_2;
				
				{
					auto const var_ref(var.ref());
					auto const var_ref_size(var_ref.size());
					auto const end(pos + var_ref_size);
					auto const var_lineno(var.lineno());
				
					// First check that there is at least one ALT that can be handled.
					if (!can_handle_variant_alts(var))
					{
						skipped_variants.insert(var_lineno);
						error_logger.log_no_supported_alts(var_lineno);
						goto loop_end_2;
					}
					
					{
						// Try to find an end position that is greater than var's position.
						auto it(end_positions.upper_bound(pos));
						auto const end_it(end_positions.cend());
						
						// If not found, add the current end position to end_positions
						// and skip the remaining checks.
						if (end_it == it)
						{
							end_positions.emplace(
								std::piecewise_construct,
								std::forward_as_tuple(end),
								std::forward_as_tuple(pos, var_lineno)
							);
							goto loop_end;
						}
						
						do
						{
							// Proper nesting since the current starting position must be greater
							// than the previous one.
							auto const other_end(it->first);
							if (end <= other_end)
								break;
							
							// Check if the potentially conflicting variant is in fact inside this one.
							auto const other_lineno(it->second.lineno);
							auto const other_pos(it->second.pos);
							if (pos == other_pos)
								goto loop_end_3;
							
							++conflict_count;
							
							// Convert starting to 1-based to get ranges like [x, y] (instead of [x, y)).
							std::cerr
							<< "Variant on line " << var_lineno << " conflicts with line " << other_lineno
							<< " ([" << 1 + pos << ", " << end << "] vs. [" << 1 + other_pos << ", " << other_end << "])." << std::endl;
							
							{
								auto const res(bad_overlaps.insert(overlap_map::value_type(other_lineno, var_lineno)));
								libbio_always_assert_msg(res.second, "Unable to insert");
							}
							
							++conflict_counts.left[other_lineno];
							++conflict_counts.left[var_lineno];
							
						loop_end_3:
							++it;
						} while (end_it != it);
					}
			
					// Add the end position.
					end_positions.emplace(
						std::piecewise_construct,
						std::forward_as_tuple(end),
						std::forward_as_tuple(pos, var_lineno)
					);
				}
				
			loop_end:
				last_position = pos;
			
			loop_end_2:
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
		
		libbio_always_assert_msg(bad_overlaps.size() == 0, "Unable to remove all conflicting variants");
		
		return conflict_count;
	}
}
