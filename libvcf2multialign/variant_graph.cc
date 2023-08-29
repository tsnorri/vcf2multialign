/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/vcf/vcf_reader.hh>
#include <vcf2multialign/transpose_matrix.hh>
#include <vcf2multialign/variant_graph.hh>

namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace vcf	= libbio::vcf;


namespace {
	
	// Boilerplate.
	// FIXME: come up with a way not to duplicate the code needed for storing field pointers.
	struct variant_format final : public vcf::variant_format
	{
		vcf::genotype_field_gt	*gt_field{};
		
		// Return a new empty instance of this class.
		virtual variant_format *new_instance() const override { return new variant_format(); }
		
		virtual void reader_did_update_format(vcf::reader &reader) override
		{
			this->assign_field_ptr("GT", gt_field);
		}
	};
	
	inline variant_format const &get_variant_format(vcf::variant const &var)
	{
		libbio_assert(var.reader()->has_assigned_variant_format());
		return static_cast <variant_format const &>(var.get_format());
	}
	
	inline variant_format const &get_variant_format(vcf::transient_variant const &var)
	{
		libbio_assert(var.reader()->has_assigned_variant_format());
		return static_cast <variant_format const &>(var.get_format());
	}
}


namespace vcf2multialign {
	
	auto variant_graph::add_node(position_type const ref_pos, position_type const aln_pos) -> node_type
	{
		reference_positions.emplace_back(ref_pos);
		aligned_positions.emplace_back(aln_pos);
		alt_edge_count_csum.emplace_back(alt_edge_count_csum.back());
		return reference_positions.size() - 1;
	}


	auto variant_graph::add_or_update_node(position_type const ref_pos, position_type const aln_pos) -> node_type
	{
		auto const last_ref_pos(reference_positions.back());
		libbio_assert_lte(last_ref_pos, ref_pos);
	
		if (last_ref_pos < ref_pos)
			return add_node(ref_pos, aln_pos);
	
		aligned_positions.back() = std::max(aligned_positions.back(), aln_pos);
		return reference_positions.size() - 1;
	}


	auto variant_graph::add_edge(std::string_view const label) -> edge_type
	{
		++alt_edge_count_csum.back();
		alt_edge_targets.emplace_back();
		alt_edge_labels.emplace_back(label);
		return alt_edge_targets.size() - 1;
	}
	
	
	void build_variant_graph(
		sequence_type const &ref_seq,
		char const *variants_path,
		char const *chr_id,
		variant_graph &graph,
		build_graph_statistics &stats,
		std::ostream *log_stream
	)
	{
		typedef variant_graph::position_type	position_type;
		typedef variant_graph::edge_type		edge_type;
		
		struct edge_destination
		{
			edge_type		edge_index{};
			position_type	position{};
		};
		
		constexpr std::size_t const path_matrix_row_col_divisor{64}; // Make sure we can transpose the matrix with the 8×8 operation.
		constexpr std::size_t const path_column_allocation{512};
		
		// FIXME: Use the BCF library when it is ready.
		// Open the variant file.
		vcf::mmap_input vcf_input;
		vcf_input.handle().open(variants_path);
		
		vcf::reader reader(vcf_input);
		
		vcf::add_reserved_info_keys(reader.info_fields());
		vcf::add_reserved_genotype_keys(reader.genotype_fields());
		
		// Read the headers.
		reader.set_variant_format(new variant_format());
		reader.read_header();
		reader.set_parsed_fields(vcf::field::ALL);
		
		graph.sample_names = reader.sample_names_by_index();
		graph.alt_edge_count_csum.emplace_back(0);
		graph.add_node(0, 0);
		
		std::uint64_t var_idx{};
		position_type aln_pos{};
		position_type prev_ref_pos{};
		bool is_first{true};
		std::vector <edge_type> edges_by_alt;
		std::vector <position_type> target_ref_positions_by_chrom_copy;
		std::vector <position_type> current_edge_targets;
		std::multimap <position_type, edge_destination> next_aligned_positions; // Aligned positions by reference position.
		
		auto add_target_nodes([&graph, &aln_pos, &prev_ref_pos, &next_aligned_positions](position_type const ref_pos){
			auto const end(next_aligned_positions.end());
			auto it(next_aligned_positions.begin());
			auto const begin(it);
			for (; it != end; ++it)
			{
				if (ref_pos < it->first)
					break;
				
				// Add the node.
				auto const dist(it->first - prev_ref_pos); // Distance from the previous node.
				aln_pos = std::max(aln_pos + dist, it->second.position);
				auto const node_idx(graph.add_or_update_node(it->first, aln_pos));
				graph.alt_edge_targets[it->second.edge_index] = node_idx;
				prev_ref_pos = it->first;
			}
			
			next_aligned_positions.erase(begin, it);
		});
		
		reader.parse(
			[
				chr_id,
				&graph,
				&stats,
				log_stream, // Pointer
				&var_idx,
				&aln_pos,
				&prev_ref_pos,
				&is_first,
				&edges_by_alt,
				&target_ref_positions_by_chrom_copy,
				&current_edge_targets,
				&next_aligned_positions,
				&add_target_nodes
			](vcf::transient_variant const &var) -> bool {
				++var_idx;
				auto const * const gt_field(get_variant_format(var).gt_field);
				
				if (var.chrom_id() != chr_id)
				{
					++stats.chr_id_mismatches;
					goto end;
				}
				
				if (!gt_field)
				{
					std::cerr << "ERROR: Variant " << var_idx << " does not have a genotype.\n";
					std::exit(EXIT_FAILURE);
				}
				
				if (is_first)
				{
					is_first = false;
					
					// Check the ploidy.
					graph.ploidy_csum.clear();
					graph.ploidy_csum.resize(1 + graph.sample_names.size(), 0);
					for (auto const &[sample_idx, sample] : rsv::enumerate(var.samples()))
					{
						auto const &gt((*gt_field)(sample));
						graph.ploidy_csum[1 + sample_idx] = graph.ploidy_csum[sample_idx] + gt.size();
					}
					
					{
						// Make sure the row count is divisible by path_matrix_row_col_divisor.
						std::size_t const path_matrix_rows(path_matrix_row_col_divisor * std::ceil(1.0 * graph.ploidy_csum.back() / path_matrix_row_col_divisor)); 
						graph.paths_by_chrom_copy_and_edge = variant_graph::path_matrix(path_matrix_rows, path_column_allocation);
					}
					
					target_ref_positions_by_chrom_copy.resize(graph.ploidy_csum.back(), 0);
				}
				
				{
					++stats.handled_variants;
					auto const ref_pos(var.zero_based_pos());
					if (! (prev_ref_pos <= ref_pos))
					{
						std::cerr << "ERROR: Variant " << var_idx << " has non-increasing position (" << prev_ref_pos << " v. " << ref_pos << ").\n";
						std::exit(EXIT_FAILURE);
					}
				
					// Add nodes if needed.
					add_target_nodes(ref_pos);
					
					// Add node for the current record.
					auto const dist(ref_pos - prev_ref_pos);
					aln_pos += dist;
					auto const node_idx(graph.add_or_update_node(ref_pos, aln_pos));
					
					// Add edges.
					// FIXME: Take END into account here?
					auto const &ref(var.ref());
					auto const &alts(var.alts());
					edges_by_alt.clear();
					edges_by_alt.resize(alts.size(), variant_graph::EDGE_MAX);
					edge_type min_edge{}, max_edge{};
					current_edge_targets.clear();
					{
						bool is_first{true};
						for (auto const &[alt_idx, alt] : rsv::enumerate(alts))
						{
							switch (alt.alt_sv_type)
							{
								case vcf::sv_type::NONE:
								case vcf::sv_type::DEL:
								{
									auto const ref_target_pos(ref_pos + ref.size());
									auto const edge_idx([&](){
										if (vcf::sv_type::NONE == alt.alt_sv_type)
										{
											auto const edge_idx(graph.add_edge(alt.alt));
											next_aligned_positions.emplace(ref_target_pos, edge_destination{edge_idx, aln_pos + alt.alt.size()});
											return edge_idx;
										}
										else
										{
											auto const edge_idx(graph.add_edge());
											next_aligned_positions.emplace(ref_target_pos, edge_destination{edge_idx, aln_pos});
											return edge_idx;
										}
									}());
									
									edges_by_alt[alt_idx] = edge_idx;
									current_edge_targets.emplace_back(ref_target_pos);
									
									if (is_first)
										min_edge = edge_idx;
									max_edge = edge_idx;
									break;
								}
								
								default:
									break;
							}
						}
					}
					
					// Check that we have enough space for the paths.
					{
						auto const ncol(graph.paths_by_chrom_copy_and_edge.number_of_columns());
						if (ncol <= max_edge)
						{
							auto const multiplier(4 + ncol / path_column_allocation);
							graph.paths_by_chrom_copy_and_edge.resize(graph.paths_by_chrom_copy_and_edge.number_of_rows() * multiplier * path_column_allocation, 0);
						}
					}
					
					// Paths.
					for (auto const &[sample_idx, sample] : rsv::enumerate(var.samples()))
					{
						auto const &gt((*gt_field)(sample));
						libbio_assert_lt(sample_idx, graph.ploidy_csum.size());
						auto const base_idx(graph.ploidy_csum[sample_idx]); // Base index for this sample.
						
						for (auto const &[chr_idx, sample_gt] : rsv::enumerate(gt))
						{
							if (0 == sample_gt.alt)
								continue;
							
							if (vcf::sample_genotype::NULL_ALLELE == sample_gt.alt)
								continue;
							
							// Check that the alternative allele was handled.
							auto const edge_idx(edges_by_alt[sample_gt.alt - 1]);
							if (variant_graph::EDGE_MAX == edge_idx)
								continue;
							
							// Check for overlapping edges for the current sample.
							auto const row_idx(base_idx + chr_idx);
							if (ref_pos < target_ref_positions_by_chrom_copy[row_idx])
							{
								if (log_stream)
									(*log_stream) << "Overlapping alternative alleles. Sample: " << graph.sample_names[sample_idx] << " chromosome copy: " << chr_idx << " current variant position: " << ref_pos << '\n';
								continue;
							}
							
							// Update the target position for this chromosome copy and the path information.
							auto const target_pos(current_edge_targets[edge_idx - min_edge]);
							target_ref_positions_by_chrom_copy[row_idx] = target_pos;
							graph.paths_by_chrom_copy_and_edge(row_idx, edge_idx) |= 1;
						}
					}
					
					prev_ref_pos = ref_pos;
				}

			end:
				if (0 == var_idx % 1'000'000)
					lb::log_time(std::cerr) << "Handled " << var_idx << " variants…\n";
				return true; // Continue parsing.
			}
		);
		
		// Add a sink node.
		{
			auto const ref_pos(ref_seq.size());
			add_target_nodes(ref_pos);
			auto const dist(ref_pos - prev_ref_pos);
			graph.add_or_update_node(ref_pos, aln_pos + dist);
		}
		
		// Remove extra entries.
		// Does not actually shrink to fit.
		{
			std::size_t const ncol(path_matrix_row_col_divisor * std::ceil(1.0 * graph.edge_count() / path_matrix_row_col_divisor));
			graph.paths_by_chrom_copy_and_edge.resize(graph.paths_by_chrom_copy_and_edge.number_of_rows() * ncol, 0);
		}
		
		graph.paths_by_chrom_copy_and_edge = transpose_matrix(graph.paths_by_chrom_copy_and_edge);
	}
}
