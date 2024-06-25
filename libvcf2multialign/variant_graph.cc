/*
 * Copyright (c) 2023-2024 Tuukka Norri
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


	struct sample_chromosome_index
	{
		std::uint32_t	sample_vcf_index{};
		std::uint32_t	sample_output_index{};
		std::uint32_t	chromosome_copy_vcf_index{};
		std::uint32_t	chromosome_copy_output_index{};

		auto to_tuple() const { return std::make_tuple(sample_vcf_index, sample_output_index, chromosome_copy_vcf_index, chromosome_copy_output_index); }
		bool operator<(sample_chromosome_index const &other) const { return to_tuple() < other.to_tuple(); }
	};
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
		build_graph_delegate &delegate
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
		
		std::string_view const ref_seq_sv{ref_seq.data(), ref_seq.size()};
		
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
		std::vector <sample_chromosome_index> included_samples;
		
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
				libbio_assert_lt(it->second.edge_index, graph.alt_edge_targets.size());
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
				&delegate,
				ref_seq_sv,
				&reader,
				&var_idx,
				&aln_pos,
				&prev_ref_pos,
				&is_first,
				&edges_by_alt,
				&target_ref_positions_by_chrom_copy,
				&current_edge_targets,
				&next_aligned_positions,
				&included_samples,
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
					
					// Check the ploidy and sample inclusion.
    				graph.ploidy_csum.clear();
					graph.ploidy_csum.resize(1 + graph.sample_names.size(), 0);
					std::vector <std::uint32_t> removed_samples; // No chromosome copies included.
					auto const &sample_names(reader.sample_names_by_index());
					std::size_t sample_idx_output{};
					for (auto const &[sample_idx_input, sample] : rsv::enumerate(var.samples()))
					{
						auto const &gt((*gt_field)(sample));
						variant_graph::ploidy_type included_count{};
						for (auto const chrom_copy_idx : rsv::iota(0U, gt.size()))
						{
							if (delegate.should_include(sample_names[sample_idx_input], chrom_copy_idx))
							{
								included_samples.emplace_back(sample_idx_input, sample_idx_output, chrom_copy_idx, included_count);
								++included_count;
							}
						}

						if (included_count)
						{
							libbio_assert_lt(1 + sample_idx_output, graph.ploidy_csum.size());
							graph.ploidy_csum[1 + sample_idx_output] = graph.ploidy_csum[sample_idx_output] + included_count;
							++sample_idx_output;
						}
						else
						{
							removed_samples.push_back(sample_idx_input);
						}
					}

					// If sample names were removed, replace the name vector.
					if (!removed_samples.empty())
					{
						libbio_assert_lte(removed_samples.size(), graph.ploidy_csum.size());
						graph.ploidy_csum.resize(graph.ploidy_csum.size() - removed_samples.size());
						
						removed_samples.push_back(UINT32_MAX);
						auto it(removed_samples.begin());
						variant_graph::label_vector new_sample_names;
						
						libbio_assert_lte((removed_samples.size() - 1), graph.sample_names.size());
						new_sample_names.reserve(graph.sample_names.size() - (removed_samples.size() - 1));
						
						for (auto &&[idx, sample_name] : rsv::enumerate(graph.sample_names))
						{
							if (idx < *it)
								new_sample_names.emplace_back(std::move(sample_name));
							else
								++it;
						}

						using std::swap;
						swap(graph.sample_names, new_sample_names);
					}
					
					{
						// Make sure the row count is divisible by path_matrix_row_col_divisor.
						std::size_t const path_matrix_rows(path_matrix_row_col_divisor * std::ceil(1.0 * graph.ploidy_csum.back() / path_matrix_row_col_divisor)); 
						graph.paths_by_edge_and_chrom_copy = variant_graph::path_matrix(path_matrix_rows, path_column_allocation);
					}
					
					libbio_assert_eq(graph.ploidy_csum.size(), 1 + graph.sample_names.size());
					libbio_assert_lte(graph.sample_names.size(), graph.ploidy_csum.back());
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
					
					// Compare to the reference.
					auto const &ref(var.ref());
					
					{
						auto const expected_ref(ref_seq_sv.substr(ref_pos, ref.size()));
						if (ref != expected_ref && !delegate.ref_column_mismatch(var_idx, var, expected_ref))
							return false;
					}
					
					// Add the edges.
					// We add even if none of the paths has the edge.
					// FIXME: Take END into account here?
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
									
									libbio_assert_lt(alt_idx, edges_by_alt.size());
									edges_by_alt[alt_idx] = edge_idx;
									current_edge_targets.emplace_back(ref_target_pos);
									
									if (is_first)
									{
										min_edge = edge_idx;
										is_first = false;
									}
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
						auto const ncol(graph.paths_by_edge_and_chrom_copy.number_of_columns());
						if (ncol <= max_edge)
						{
							auto const multiplier(4 + ncol / path_column_allocation);
							graph.paths_by_edge_and_chrom_copy.resize(graph.paths_by_edge_and_chrom_copy.number_of_rows() * multiplier * path_column_allocation, 0);
						}
					}
					
					// Paths.
					for (auto const &sample_chr_idx : included_samples)
					{
						auto const sample_idx_input(sample_chr_idx.sample_vcf_index);
						auto const sample_idx_output(sample_chr_idx.sample_output_index);
						auto const chr_idx_input(sample_chr_idx.chromosome_copy_vcf_index);
						auto const chr_idx_output(sample_chr_idx.chromosome_copy_output_index);

						auto const &sample(var.samples()[sample_idx_input]);
						auto const &gt((*gt_field)(sample));
						libbio_assert_lt(sample_idx_output, graph.ploidy_csum.size());
						auto const base_idx(graph.ploidy_csum[sample_idx_output]); // Base index for this sample.
						libbio_assert_lt(chr_idx_input, gt.size());
						auto const &sample_gt(gt[chr_idx_input]);

						if (0 == sample_gt.alt)
							continue;
						
						if (vcf::sample_genotype::NULL_ALLELE == sample_gt.alt)
							continue;
						
						// Check that the alternative allele was handled.
						libbio_assert_lt(sample_gt.alt - 1, edges_by_alt.size());
						auto const edge_idx(edges_by_alt[sample_gt.alt - 1]);
						if (variant_graph::EDGE_MAX == edge_idx)
							continue;
						
						// Check for overlapping edges for the current sample.
						auto const row_idx(base_idx + chr_idx_output);
						libbio_assert_lt(row_idx, target_ref_positions_by_chrom_copy.size());
						if (ref_pos < target_ref_positions_by_chrom_copy[row_idx])
						{
							delegate.report_overlapping_alternative(
								var.lineno() + reader.last_header_lineno(),
								ref_pos,
								var.id(),
								reader.sample_names_by_index()[sample_idx_input],
								chr_idx_input,
								sample_gt.alt
							);
						}
							
						// Update the target position for this chromosome copy and the path information.
						libbio_assert_lt(edge_idx - min_edge, current_edge_targets.size());
						auto const target_pos(current_edge_targets[edge_idx - min_edge]);
						target_ref_positions_by_chrom_copy[row_idx] = target_pos;
						graph.paths_by_edge_and_chrom_copy(row_idx, edge_idx) |= 1;
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
			graph.paths_by_edge_and_chrom_copy.resize(graph.paths_by_edge_and_chrom_copy.number_of_rows() * ncol, 0);
		}
		
		graph.paths_by_chrom_copy_and_edge = transpose_matrix(graph.paths_by_edge_and_chrom_copy);
	}
}
