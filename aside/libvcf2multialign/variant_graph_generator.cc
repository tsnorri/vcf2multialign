/*
 * Copyright (c) 2019-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/algorithm.hh>
#include <libbio/bits.hh>
#include <libbio/vcf/variant_printer.hh>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/transform.hpp>
#include <vcf2multialign/variant_graph/variant_graph_generator.hh>
#include <vcf2multialign/utility/can_handle_variant_alts.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace vcf	= libbio::vcf;


namespace {

	template <typename t_container>
	class back_emplace_iterator
	{
	public:
		typedef t_container	container_type;
		typedef ptrdiff_t	difference_type;
		
	protected:
		container_type	*m_container{};
		
	public:
		back_emplace_iterator() = default;
		
		explicit back_emplace_iterator(t_container &container):
			m_container(&container)
		{
		}
		
		template <typename t_arg>
		back_emplace_iterator &operator=(t_arg &&arg)
		{
			m_container->emplace_back(std::forward <t_arg>(arg));
			return *this;
		}
		
		back_emplace_iterator &operator*() { return *this; }
		back_emplace_iterator &operator++() { return *this; }
		back_emplace_iterator &operator++(int) { return *this; }
	};
	
	
	template <typename t_container>
	inline back_emplace_iterator <t_container> back_emplacer(t_container &container)
	{
		return back_emplace_iterator(container);
	}
}


namespace vcf2multialign { namespace variant_graphs {
	
	void variant_graph_generator::update_sample_names()
	{
		m_sample_names = this->vcf_reader().sample_names_by_index();
	}
	
	
	void variant_graph_generator::generate_graph_setup()
	{
		auto &reader(this->vcf_reader());
		
		// Update the delegate in case *this has been moved.
		m_sample_sorter.set_delegate(*this);
		
		update_sample_names();
		m_graph.clear();
		
		m_processed_count.store(0, std::memory_order_relaxed);
		reader.set_parsed_fields(vcf::field::ALL);
	}
	
	
	void variant_graph_single_pass_generator::generate_graph(
		std::vector <std::string> const &field_names_for_filter_by_assigned,
		bool const should_start_from_current_variant
	)
	{
		auto &reader(this->vcf_reader());
		auto &delegate(this->delegate());
		
		generate_graph_setup();
		
		// Determine the fields used for filtering.
		std::vector <vcf::info_field_base *> filter_by_assigned;
		this->fill_filter_by_assigned(field_names_for_filter_by_assigned, filter_by_assigned, delegate);
		
		std::size_t prev_overlap_end_pos(0);
		std::size_t overlap_end_pos(0);
		bool is_first{true};
		
		// Process the variants.
		auto handling_callback(
			[
				this,
				&delegate,
				&filter_by_assigned,
				&overlap_end_pos,
				&prev_overlap_end_pos,
				&is_first
			]
			(vcf::transient_variant const &var) -> bool
			{
				auto const lineno(var.lineno());
				auto const var_pos(var.zero_based_pos());
				auto const var_end(vcf::variant_end_pos(var, *m_end_field));
				libbio_always_assert_lte(var_pos, var_end);
				
				switch (this->check_variant(var, filter_by_assigned, delegate))
				{
					case variant_check_status::PASS:
						break;
					case variant_check_status::ERROR:
						goto end;
					case variant_check_status::FATAL_ERROR:
						return false;
				}
				
				// check_variant() filters by chromosome name, so we only need to prepare the sample sample indexer once.
				if (is_first)
				{
					is_first = false;
					
					libbio_assert(get_variant_format(var).gt);
					m_end_positions_by_sample.clear();
					m_end_positions_by_sample.resize(m_sample_indexer.total_samples(), 0);
				}
				
				// Check for a suitable subgraph starting position.
				// process_subgraph() actually does a similar check in order to
				// handle overlaps within a subgraph.
				if (overlap_end_pos + m_minimum_bridge_length <= var_pos)
				{
					// Handle the subgraph.
					process_subgraph(prev_overlap_end_pos);
					prev_overlap_end_pos = overlap_end_pos;
				}
				
				// Add to the stack.
				overlap_end_pos = std::max(overlap_end_pos, var_end);
				m_subgraph_variants.emplace_back(var);
				delegate.variant_processor_found_matching_variant(var);
				
			end:
				++m_processed_count;
				return true;
			}
		);
		
		bool should_continue(true);
		if (should_start_from_current_variant)
			should_continue = handling_callback(reader.current_variant());
		if (should_continue)
			reader.parse(handling_callback);
		
		// Process the remaining variants.
		process_subgraph(prev_overlap_end_pos);
		
		// Copy the sample names and finalize.
		m_graph.sample_names() = m_sample_names;
		finalize_graph();
	}
	
	
	void variant_graph_precalculated_generator::generate_graph(bool const should_start_from_current_variant)
	{
		auto &reader(this->vcf_reader());
		auto &delegate(this->delegate());
		
		generate_graph_setup();
		
		// Prepare using the preprocessing result.
		auto const &preprocessing_res(*m_preprocessing_result);
		auto const cut_pos_rng(preprocessing_res.positions | rsv::tail);
		auto cut_pos_it(ranges::begin(cut_pos_rng));
		auto const cut_pos_end(ranges::end(cut_pos_rng));
		auto line_no_it(ranges::begin(preprocessing_res.handled_line_numbers));
		auto const line_no_end(ranges::end(preprocessing_res.handled_line_numbers));
		std::size_t overlap_end_pos(0);
		std::size_t prev_overlap_end_pos(0);
		libbio_always_assert_neq(cut_pos_it, cut_pos_end);
		bool is_first{true};
		
		// Process the variants.
		auto handling_callback(
			[
				this,
				&cut_pos_it,
				cut_pos_end,
				&line_no_it,
				line_no_end,
				&overlap_end_pos,
				&prev_overlap_end_pos,
				&is_first
			](vcf::transient_variant const &var) -> bool
			{
				auto const lineno(var.lineno());
				libbio_always_assert_lte(lineno, *line_no_it);
				if (lineno == *line_no_it)
				{
					auto const var_pos(var.zero_based_pos());
					auto const var_end(vcf::variant_end_pos(var, *m_end_field));
					libbio_always_assert_lte(var_pos, var_end);
					libbio_always_assert_lte(var_pos, *cut_pos_it);
					
					if (is_first)
					{
						is_first = false;
						libbio_assert(get_variant_format(var).gt);
						m_end_positions_by_sample.clear();
						m_end_positions_by_sample.resize(m_sample_indexer.total_samples(), 0);
					}
					
					if (var_pos == *cut_pos_it)
					{
						// Handle the subgraph.
						process_subgraph(prev_overlap_end_pos);
						prev_overlap_end_pos = overlap_end_pos;
						
						++cut_pos_it;
						libbio_always_assert_neq(cut_pos_it, cut_pos_end);
					}
					
					// Add to the stack.
					overlap_end_pos = std::max(overlap_end_pos, var_end);
					m_subgraph_variants.emplace_back(var);
					
					++m_processed_count;
					++line_no_it;
					
					auto const should_continue(line_no_it != line_no_end);
					return should_continue;
				}
				return true;
			}
		);
		
		bool should_continue(true);
		if (should_start_from_current_variant)
			should_continue = handling_callback(reader.current_variant());
		if (should_continue)
			reader.parse(handling_callback);
		
		// Process the remaining variants.
		process_subgraph(prev_overlap_end_pos);
		
		// Copy the sample names and finalize.
		m_graph.sample_names() = m_sample_names;
		finalize_graph();
	}
	
	
	void variant_graph_generator::combine_subgraph_variants_by_pos_and_ref()
	{
		// Check consecutive variants in m_subgraph_variants for matching POS and REF values.
		// If any are found, set the GT values in the left hand side and remove the right hand side variant.
		// (No attempt is made to merge the other INFO and FORMAT values.)

		if (m_subgraph_variants.size() < 2)
			return;

		auto it(m_subgraph_variants.begin() + 1);
		while (it != m_subgraph_variants.end())
		{
			auto &prev_var(*(it - 1));
			auto const &curr_var(*it);

			if (prev_var.zero_based_pos() == curr_var.zero_based_pos())
			{
				if (prev_var.ref() == curr_var.ref())
				{
					// Move the ALT and GT values.
					auto &prev_alts(prev_var.alts());
					auto const *prev_gt_field(get_variant_format(prev_var).gt);
					auto const *curr_gt_field(get_variant_format(curr_var).gt);
					for (auto const &[i, curr_alt] : rsv::enumerate(curr_var.alts()))
					{
						auto const curr_alt_idx(1 + i);
						auto const prev_alt_it(std::find(prev_alts.begin(), prev_alts.end(), curr_alt));
						std::size_t new_alt_idx(0);
						if (prev_alts.end() == prev_alt_it)
						{
							// ALT not found in the previous variant.
							// Add the ALT to the end.
							prev_alts.push_back(curr_alt);
							new_alt_idx = prev_alts.size();
						}
						else
						{
							// ALT found in the previous variant.
							new_alt_idx = 1 + std::distance(prev_alts.begin(), prev_alt_it);
						}
						// Copy the GT values.
						// We assume that the ploidy does not change.
						for (auto &&[prev_sample, curr_sample] : rsv::zip(prev_var.samples(), curr_var.samples()))
						{
							auto &prev_gt_vec((*prev_gt_field)(prev_sample));
							auto const curr_gt_vec((*curr_gt_field)(curr_sample));
							for (auto &&[prev_gt_val, curr_gt_val] : rsv::zip(prev_gt_vec, curr_gt_vec))
							{
								if (curr_gt_val.alt == curr_alt_idx)
								{
									// GT value matches the ALT index, copy.
									prev_gt_val.alt = new_alt_idx;
								}
							}
						}
					}

					// Remove the variant and continue.
					it = m_subgraph_variants.erase(it);
					continue;
				}
			}

			++it;
		}
	}
	
	
	void variant_graph_generator::reset_genotype_values_for_samples_with_overlapping_variants()
	{
		// Check each sample for overlapping variants and set GT to zero if necessary.
		// Also zero ALT indices that cannot be handled, as well as null alleles.
		
		auto const total_samples(m_sample_indexer.total_samples());
		libbio_assert_lt(0, total_samples);
		for (auto &var : m_subgraph_variants)
		{
			auto const pos(var.zero_based_pos());
			auto const end_pos(vcf::variant_end_pos(var, *m_end_field));
			auto const &alts(var.alts());
			auto &samples(var.samples());
			auto const *gt_field(get_variant_format(var).gt);
			
			// Determine which ALTs can be handled.
			m_can_handle_alt.clear();
			m_can_handle_alt.resize(alts.size(), 0);
			for (auto const &[i, alt] : rsv::enumerate(var.alts()))
				m_can_handle_alt[i] = can_handle_variant_alt(alt);
			
			// Check the samples.
			for (auto const sample_idx : rsv::iota(0U, total_samples))
			{
				libbio_assert_lt_msg(sample_idx, m_end_positions_by_sample.size(), "sample_idx: ", sample_idx, " m_end_positions_by_sample.size(): ", m_end_positions_by_sample.size());
				auto const [donor_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(sample_idx);
				
				auto &sample(samples[donor_idx]);
				auto &gt((*gt_field)(sample));
				auto const alt_idx(gt[chr_idx].alt); // bit-field, cannot be referenced.
				
				if (vcf::sample_genotype::NULL_ALLELE == alt_idx)
					gt[chr_idx].alt = 0;
				else if (alt_idx)
				{
					if (m_can_handle_alt[alt_idx - 1])
					{
						// ALT can be handled. Check the position.
						if (m_end_positions_by_sample[sample_idx] <= pos)
							m_end_positions_by_sample[sample_idx] = end_pos;
						else
						{
							// Sample has overlapping variants.
							// For 1000G variants setting GT to zero is ok b.c. the samples
							// do not have variant probabilities.
							gt[chr_idx].alt = 0;
						}
					}
					else
					{
						// ALT cannot be handled, set to zero.
						gt[chr_idx].alt = 0;
					}
				}
			}
		}
	}
	
	
	std::tuple <std::size_t, std::size_t> variant_graph_generator::create_subgraph_and_nodes()
	{
		auto const nodes_with_alts(ranges::count_if(m_sorted_nodes, [](auto const &node){ return 0 < node.alt_edge_count; }));
		
		auto const sample_count(m_sample_indexer.total_samples());
		libbio_assert_lt(0, sample_count);
		auto const path_count(m_sample_sorter.path_count());
		this->delegate().variant_graph_generator_will_handle_subgraph(m_subgraph_variants.front(), m_subgraph_variants.size(), path_count);
		
		auto &first_node(m_sorted_nodes.front());
		auto const &[first_node_idx, first_alt_edge_start] = m_graph.add_main_node_if_needed(first_node.ref_position, first_node.alt_edge_count);
		first_node.node_index = first_node_idx;
		first_node.alt_edge_start = first_alt_edge_start;
	
		auto const subgraph_start_pos(first_node.ref_position);
		auto const subgraph_idx(m_graph.add_subgraph(first_node_idx, sample_count, path_count));
		auto max_alt_edge_count(first_node.alt_edge_count);
		
		for (auto &node : m_sorted_nodes | rsv::tail)
		{
			auto const &[node_idx, alt_edge_start] = m_graph.add_main_node(node.ref_position, node.alt_edge_count);
			node.node_index = node_idx;
			node.alt_edge_start = alt_edge_start;
			max_alt_edge_count = std::max(max_alt_edge_count, node.alt_edge_count);
		}
	
		m_graph.setup_path_edges_for_current_subgraph(max_alt_edge_count, nodes_with_alts, path_count);
	
		return {subgraph_idx, first_node_idx};
	}
	
	
	// Retrieve the accumulated group of variants and pass them to the worker thread for processing.
	void variant_graph_generator::process_subgraph(std::size_t const prev_overlap_end_pos)
	{
		// Fast path: empty stack.
		if (m_subgraph_variants.empty())
			return;
	
		auto const &ref_positions(m_graph.ref_positions());
		auto const &alt_edge_count_csum(m_graph.alt_edge_count_csum());
		auto const &alt_edge_targets(m_graph.alt_edge_targets());
		auto const &alt_edge_labels(m_graph.alt_edge_labels());
	
		// Slow path: possibly overlapping variants.
		// Combine variants that have matching POS and REF.
		combine_subgraph_variants_by_pos_and_ref();
		
		// Fix bad overlaps, i.e. for overlapping variants a and b, sample s has both a and b.
		reset_genotype_values_for_samples_with_overlapping_variants();
	
		m_sample_sorter.prepare_for_next_subgraph();
	
		// Sort by variant and ALT.
		for (auto const &var : m_subgraph_variants)
		{
			for (std::size_t i(0), count(var.alts().size()); i < count; ++i)
				m_sample_sorter.sort_by_variant_and_alt(var, 1 + i);
	
			// FIXME: consider handling null alleles.
		}
	
		// Determine the nodes to be created.
		// Begin by listing all the distinct start and end positions.
		m_sorted_nodes.clear();
		m_start_positions.clear();
		m_end_positions.clear();
		for (auto const &var : m_subgraph_variants)
		{
			auto const ref_pos(var.zero_based_pos());
			auto const end_pos(vcf::variant_end_pos(var, *m_end_field));
			m_start_positions.push_back(ref_pos);
			m_end_positions.push_back(end_pos);
		}

		{
			auto const end(m_start_positions.end());
			auto const it(std::is_sorted_until(m_start_positions.begin(), end));
			if (it != end)
			{
				std::cerr << "ERROR: The following variants are not sorted by position in the input.\n";
				for (auto const &var : m_subgraph_variants)
					vcf::output_vcf(std::cerr, var);
				std::cerr << "First position not in sorted order:" << (1 + *it) << "\n";
				// Throw instead of std::exit?
				std::exit(EXIT_FAILURE);
			}
		}

		ranges::sort(m_end_positions);
	
		// Get unique positions and merge.
		ranges::set_union(m_start_positions | rsv::unique, m_end_positions | rsv::unique, back_emplacer(m_sorted_nodes));
	
		// Store the ALT counts.
		for (auto const &var : m_subgraph_variants)
		{
			auto const ref_pos(var.zero_based_pos());
			auto const handled_alt_count(ranges::count_if(var.alts(), [](auto const &alt){ return can_handle_variant_alt(alt); }));
			auto const it(ranges::partition_point(m_sorted_nodes, [ref_pos](detail::generator_node_description const &node){
				return node.ref_position < ref_pos;
			}));
	
			libbio_assert_neq(it, m_sorted_nodes.end());
			it->alt_edge_count += handled_alt_count;
		}
	
		// Create the subgraph and the nodes.
		auto const [subgraph_idx, first_node_idx] = create_subgraph_and_nodes();
	
		// Path number for each sample.
		{
			auto const &paths_by_sample(m_sample_sorter.paths_by_sample());
			auto &paths_by_subgraph(m_graph.sample_paths());
			auto &current_subgraph_paths(paths_by_subgraph[subgraph_idx]);
			for (auto &&[src, dst] : rsv::zip(paths_by_sample, current_subgraph_paths))
				dst |= src;
		}
	
		// Connect the ALT edges.
		// Also get a representative sample for each path and store its fixed ALT indices over the variants with the corresponding path.
		{
			libbio_assert(!m_sorted_nodes.empty());
			auto &alt_edge_targets(m_graph.alt_edge_targets());
			auto &alt_edge_labels(m_graph.alt_edge_labels());
			auto node_it(m_sorted_nodes.begin());
			auto const node_end(m_sorted_nodes.end());
			std::size_t var_idx_distinct_pos(0); // Index variants that have handled ALTs within a subgraph. (Using filter and enumerate on m_subgraph_variants would be an alternative.) Count only distinct positions.
			std::size_t node_alt_edge_start(m_sorted_nodes.front().alt_edge_start); // Index of the first ALT edge for the current node.
			auto var_alt_edge_start(node_alt_edge_start); // Index of the first ALT edge for the current pair of variant and node.
			auto alt_edge_index(node_alt_edge_start); // Index of the current ALT edge.
		
			// Determine the representatives by path.
			path_sorted_variant psv;
			psv.set_paths_by_sample(m_sample_sorter.paths_by_sample());
			psv.reserve_memory_for_representatives(m_sample_sorter.path_count());
			psv.determine_representatives_for_each_sample();
			auto &dst_path_edges(m_graph.path_edges()[subgraph_idx]);
			
			auto prev_ref_pos(m_subgraph_variants.front().zero_based_pos());
			for (auto const &var : m_subgraph_variants)
			{
				auto const ref_pos(var.zero_based_pos());
				auto const end_pos(vcf::variant_end_pos(var, *m_end_field));
			
				// Prepare a cumulative sum of unhandled ALTs.
				m_unhandled_alt_csum.clear();
				auto const &var_alts(var.alts());
				auto const alt_count(var_alts.size());
				m_unhandled_alt_csum.resize(1 + alt_count);
				for (std::size_t i(0); i < alt_count; ++i)
				{
					auto const &alt(var_alts[i]);
					bool const can_handle(can_handle_variant_alt(alt));
					m_unhandled_alt_csum[1 + i] = m_unhandled_alt_csum[i] + (can_handle ? 0 : 1);
				}
				
				// Check if the current variant has any ALTs that could be handled.
				if (m_unhandled_alt_csum.back() == alt_count)
					continue;
				
				// Update the variant index.
				if (ref_pos != prev_ref_pos)
					++var_idx_distinct_pos;
				
				// Find the node description for the position of the current variant.
				if (node_it->ref_position != ref_pos)
				{
					do
					{
						++node_it;
						libbio_assert_neq(node_it, node_end);
					}
					while (node_it->ref_position != ref_pos);
				
					// If the node was changed, reset the ALT edge counters.
					node_alt_edge_start = node_it->alt_edge_start;
					var_alt_edge_start = node_alt_edge_start;
					alt_edge_index = node_alt_edge_start;
				}
			
				auto const node_idx(node_it->node_index);
			
				// Find the destination node index.
				auto const dst_node_it(ranges::partition_point(ref_positions, [end_pos](auto const &pos){
					return pos < end_pos;
				}));
				libbio_assert_neq(dst_node_it, ref_positions.end());
				auto const dst_node_idx(std::distance(ref_positions.begin(), dst_node_it) - 1); // ref_positions has 1-based indexing.
			
				// Connect the ALT edges.
				{
					std::size_t max_alt_len(0);
					for (auto const &alt : var.alts())
					{
						if (can_handle_variant_alt(alt))
						{
							libbio_assert_lte(dst_node_idx, ref_positions.size() - 2);
							alt_edge_targets[alt_edge_index] = dst_node_idx;
							switch (alt.alt_sv_type)
							{
								// Handled ALTs:
								case vcf::sv_type::NONE:
									alt_edge_labels[alt_edge_index] = alt.alt;
									break;
								case vcf::sv_type::UNKNOWN:
								case vcf::sv_type::DEL:
								case vcf::sv_type::DEL_ME:
									alt_edge_labels[alt_edge_index] = "";
									break;
								default:
									libbio_fail("Unexpected ALT type.");
									break;
							}
							++alt_edge_index;
						}
					}
				}
			
				// Set the path edges.
				auto const *gt_field(get_variant_format(var).gt);
				auto const &representatives(psv.representatives_by_path());
				for (auto const &[path_idx, sample_idx] : rsv::enumerate(representatives))
				{
					libbio_assert_neq(SAMPLE_NUMBER_MAX, sample_idx);
					auto const [donor_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(sample_idx);
					auto const &sample(var.samples()[donor_idx]);
					auto const &gt((*gt_field)(sample));
					auto const alt_idx(gt[chr_idx].alt);
				
					if (alt_idx && vcf::sample_genotype::NULL_ALLELE != alt_idx)
					{
						libbio_assert_lt(alt_idx - 1, var.alts().size());
						libbio_assert(can_handle_variant_alt(var.alts()[alt_idx - 1]));
					
						auto const fixed_alt_idx(var_alt_edge_start - node_alt_edge_start + alt_idx - m_unhandled_alt_csum[alt_idx - 1]);
						libbio_assert_lte_msg(
							fixed_alt_idx, alt_edge_count_csum[1 + node_idx] - alt_edge_count_csum[node_idx],
							"ALT index ", fixed_alt_idx, " refers to a value greater than the number of ALTs, ", 
							(alt_edge_count_csum[1 + node_idx] - alt_edge_count_csum[node_idx]),
							". Variant line number: ", var.lineno(), ", original ALT index: ", alt_idx
						);
						libbio_assert_eq(0, dst_path_edges(var_idx_distinct_pos, path_idx));
						dst_path_edges(var_idx_distinct_pos, path_idx) |= fixed_alt_idx;
					}
				}
			
				// Update the counters.
				auto const handled_alt_count(alt_count - m_unhandled_alt_csum.back());
				var_alt_edge_start += handled_alt_count;
			
				prev_ref_pos = ref_pos;
				
				libbio_assert_eq(alt_edge_index, var_alt_edge_start);
			}
		}
		
		m_graph.update_aligned_ref_positions_from_node(first_node_idx);
		
		m_subgraph_variants.clear();
	}	
	
	
	void variant_graph_generator::finalize_graph()
	{
		auto &ref_positions(m_graph.ref_positions());
		auto const ref_size(reference().size());
		if (ref_positions.back() < ref_size)
		{
			auto &aln_positions(m_graph.aligned_ref_positions());
			auto const [node_idx, alt_edge_start_idx] = m_graph.add_main_node(ref_size, 0);
			
			libbio_assert_lte(2, ref_positions.size());
			auto const idx(ref_positions.size() - 2);
			auto const prev_ref_pos(ref_positions[idx]);
			auto const prev_aln_pos(aln_positions[idx]);
			auto const length(ref_size - prev_ref_pos);
			aln_positions.back() = prev_aln_pos + length;
		}
	}
}}
