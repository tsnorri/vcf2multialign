/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/find_subgraph_starting_points.hh>
#include <vcf2multialign/util.hh>


namespace vcf2multialign {

	void find_subgraph_starting_points(
		vcf_reader &reader,
		variant_set const &skipped_variants,
		variant_set /* out */ &subgraph_starting_points
	)
	{
		std::size_t last_position(0);
		std::size_t current_subgraph_end(0);
		
		reader.reset();
		reader.set_parsed_fields(vcf_field::ALT);
		bool should_continue(false);
		do {
			reader.fill_buffer();
			should_continue = reader.parse(
				[
					&skipped_variants,
					&subgraph_starting_points,
					&last_position,
					&current_subgraph_end
				]
				(transient_variant const &var)
				-> bool
			{
				// Verify that the positions are in increasing order.
				auto const pos(var.zero_based_pos());

				always_assert(last_position <= pos, "Positions not in increasing order");
				
				auto const var_lineno(var.lineno());
				auto const var_ref(var.ref());
				auto const var_ref_size(var_ref.size());
				auto const end(pos + var_ref_size);
				
				// First check that there is at least one variant that can be handled.
				if (skipped_variants.count(var_lineno))
					goto loop_end_2;
				
				// If we've moved past the previous subgraph, record the current starting point.
				if (current_subgraph_end <= var_lineno)
					subgraph_starting_points.insert(var_lineno);
				
				// If the current variant extends the current subgraph, store the new ending point.
				// FIXME: this shouldn't happen with our current nesting rules.
				if (current_subgraph_end < end)
					current_subgraph_end = end;
				
			loop_end:
				last_position = pos;
			
			loop_end_2:
				return true;
			});
		} while (should_continue);
	}
}
