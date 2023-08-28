/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/int_matrix.hh>
#include <libbio/subprocess.hh>
#include <libbio/vcf/vcf_reader.hh>
#include <map>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <string>
#include <string_view>
#include <vcf2multialign/transpose_matrix.hh>
#include <vector>
#include "cmdline.h"

namespace ios	= boost::iostreams;
namespace lb	= libbio;
namespace rsv	= ranges::views;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace {
	
	typedef std::vector <char>	sequence_vector;
	
	
	struct variant_graph
	{
		typedef std::uint64_t				position_type;	// FIXME: is std::uint32_t enough?
		typedef std::uint64_t				node_type;
		typedef std::uint64_t				edge_type;
		typedef std::uint32_t				sample_type;
		typedef std::uint32_t				ploidy_type;
		typedef std::string					label_type;
		typedef std::vector <position_type>	position_vector;
		typedef std::vector <node_type>		node_vector;
		typedef std::vector <edge_type>		edge_vector;
		typedef std::vector <label_type>	label_vector;
		typedef std::vector <ploidy_type>	ploidy_csum_vector;
		typedef lb::bit_matrix				path_matrix;
		
		constexpr static inline edge_type const EDGE_MAX{std::numeric_limits <edge_type>::max()};
		constexpr static inline sample_type const SAMPLE_MAX{std::numeric_limits <sample_type>::max()};
		
		position_vector						reference_positions;			// Reference positions by node number.
		position_vector						aligned_positions;				// MSA positions by node number.
		node_vector							alt_edge_targets;				// ALT edge targets by edge number.
		edge_vector							alt_edge_count_csum;			// Cumulative sum of ALT edge counts by 1-based node number.
		label_vector						alt_edge_labels;				// ALT edge labels by edge number.
		path_matrix							paths_by_chrom_copy_and_edge;	// Edges on rows, chromosome copies (samples multiplied by ploidy) in columns. (Vice-versa when constructing.)
		
		label_vector						sample_names;					// Sample names by sample index. FIXME: In case we have variant_graph ->> chromosome at some point, this should be in the graph.
		ploidy_csum_vector					ploidy_csum;					// Cumulative sum of ploidies by 1-based sample number (for this chromosome).
		
		node_type node_count() const { return reference_positions.size(); }
		edge_type edge_count() const { return alt_edge_targets.size(); }
		
		std::pair <edge_type, edge_type> edge_range_for_node(node_type const &node_idx) const { return {alt_edge_count_csum[node_idx], alt_edge_count_csum[1 + node_idx]}; }
		
		ploidy_type sample_ploidy(sample_type const sample_idx) const { return ploidy_csum[1 + sample_idx] - ploidy_csum[sample_idx]; }
		
		node_type add_node(position_type const ref_pos, position_type const aln_pos);
		node_type add_or_update_node(position_type const ref_pos, position_type const aln_pos);
		edge_type add_edge(std::string_view const label = std::string_view{});
	};
	
	
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
	
	
	struct build_graph_statistics
	{
		std::uint64_t	handled_variants{};
		std::uint64_t	chr_id_mismatches{};
	};
	
	
	void build_variant_graph(sequence_vector const &ref_seq, char const *variants_path, char const *chr_id, variant_graph &graph, build_graph_statistics &stats)
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
					return true; // Continue parsing.
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
										next_aligned_positions.emplace(ref_target_pos, edge_destination{edge_idx, ref_pos + alt.alt.size()});
										return edge_idx;
									}
									else
									{
										auto const edge_idx(graph.add_edge());
										next_aligned_positions.emplace(ref_target_pos, edge_destination{edge_idx, ref_pos});
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
						auto const multiplier(1 + ncol / path_column_allocation);
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
							std::cout << "Overlapping alternative alleles. Sample: " << graph.sample_names[sample_idx] << " chromosome copy: " << chr_idx << " current variant position: " << ref_pos << '\n';
							continue;
						}
						
						// Update the target position for this chromosome copy and the path information.
						auto const target_pos(current_edge_targets[edge_idx - min_edge]);
						target_ref_positions_by_chrom_copy[row_idx] = target_pos;
						graph.paths_by_chrom_copy_and_edge(row_idx, edge_idx) |= 1;
					}
				}
				
				prev_ref_pos = ref_pos;
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
		
		v2m::transpose_matrix(graph.paths_by_chrom_copy_and_edge);
	}
	
	
	void output_graphviz(
		sequence_vector const &ref_seq_,
		variant_graph const &graph,
		lb::file_handle &fh
	)
	{
		std::string_view const ref_seq(ref_seq_.data(), ref_seq_.size());
		
		lb::file_ostream stream;
		stream.open(fh.get(), ios::never_close_handle);
		stream.exceptions(std::ostream::badbit);
		
		stream << "digraph variants {\n";
		stream << "\trankdir = LR;\n";
		stream << "\trank = same;\n";
		
		typedef variant_graph::position_type	position_type;
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::edge_type		edge_type;
		
		// Nodes.
		for (auto const &[node, ref_pos, aln_pos] : rsv::zip(rsv::iota(0), graph.reference_positions, graph.aligned_positions))
			stream << '\t' << node << " [shape = Mrecord, label = \"" << node << " | " << ref_pos << " | " << aln_pos << "\"];\n";
		stream << '\n';
		
		// REF edges.
		for (auto const &[node, range] : rsv::enumerate(graph.reference_positions | rsv::sliding(2)))
		{
			auto const lb(range[0]);
			auto const rb(range[1]);
			auto const label(ref_seq.substr(lb, rb - lb));

			// FIXME: Handle special characters in the label?
			stream << '\t' << node << " -> " << (node + 1) << " [label = \"";
			if (label.size() <= label.size())
				stream << label;
			else
				stream << label.substr(0, 10) << "…" << label.substr(label.size() - 10, 10) << " (" << label.size() << ')';
			stream << "\", penwidth = 2.0];\n";
		}
		stream << '\n';
		
		// ALT edges.
		for (auto const &[src_node, edge_range] : rsv::enumerate(graph.alt_edge_count_csum | rsv::sliding(2)))
		{
			auto const edge_lb(edge_range[0]);
			auto const edge_rb(edge_range[1]);
			
			for (edge_type edge_idx(edge_lb); edge_idx < edge_rb; ++edge_idx)
				stream << '\t' << src_node << " -> " << graph.alt_edge_targets[edge_idx] << " [label = \"" << graph.alt_edge_labels[edge_idx] << "\"];\n";
		}
		stream << "}\n";
	}
	
	
	void output_sequence_(
		sequence_vector const &ref_seq,
		variant_graph const &graph,
		variant_graph::sample_type const sample_idx,
		variant_graph::ploidy_type const chr_copy_idx,
		lb::file_handle &fh
	)
	{
		lb::file_ostream stream;
		stream.open(fh.get(), ios::never_close_handle);
		stream.exceptions(std::ostream::badbit);
		
		typedef variant_graph::position_type	position_type;
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::edge_type		edge_type;
		
		position_type ref_pos{};
		position_type aln_pos{};
		position_type next_ref_pos{};
		position_type next_aln_pos{};
		node_type current_node{};
		auto const limit(graph.node_count() - 1);
		auto const chr_copy_idx_(variant_graph::SAMPLE_MAX == sample_idx ? 0 : graph.ploidy_csum[sample_idx] + chr_copy_idx);
		while (current_node < limit)
		{
			std::size_t label_size{};
			if (variant_graph::SAMPLE_MAX != sample_idx) // Always follow REF edges if outputting the aligned reference.
			{
				auto const &[edge_lb, edge_rb] = graph.edge_range_for_node(current_node);
				for (edge_type edge_idx(edge_lb); edge_idx < edge_rb; ++edge_idx)
				{
					if (graph.paths_by_chrom_copy_and_edge(edge_idx, chr_copy_idx_))
					{
						// Found an ALT edge to follow.
						auto const target_node(graph.alt_edge_targets[edge_idx]);
						auto const &label(graph.alt_edge_labels[edge_idx]);
						next_ref_pos = graph.reference_positions[target_node];
						next_aln_pos = graph.aligned_positions[target_node];
						libbio_assert_lte(label.size(), next_aln_pos - aln_pos);
						stream << label;
						current_node = target_node;
						label_size = label.size();
						goto continue_loop;
					}
				}
			}
			
			{
				next_ref_pos = graph.reference_positions[current_node + 1];
				next_aln_pos = graph.aligned_positions[current_node + 1];
				std::string_view const ref_part(ref_seq.data() + ref_pos, next_ref_pos - ref_pos);
				stream << ref_part;
				label_size = ref_part.size();
				++current_node;
			}
			
		continue_loop:
			std::fill_n(std::ostreambuf_iterator <char>(stream), next_aln_pos - aln_pos - label_size, '-');
			ref_pos = next_ref_pos;
			aln_pos = next_aln_pos;
		}
	}
	
	
	void output_sequence(
		sequence_vector const &ref_seq,
		variant_graph const &graph,
		variant_graph::sample_type const sample_idx,
		variant_graph::ploidy_type const chr_copy_idx,
		char const * const pipe_cmd,
		char const * const dst_name
	)
	{
		if (pipe_cmd)
		{
			auto proc(lb::subprocess <lb::subprocess_handle_spec::STDIN>::subprocess_with_arguments({pipe_cmd, dst_name}));
			auto &fh(proc.stdin_handle());
			
			output_sequence_(ref_seq, graph, sample_idx, chr_copy_idx, fh);
			
			auto const res(proc.close());
			auto const close_status(std::get <0>(res));
			auto const exit_status(std::get <1>(res));
			if (! (lb::process_handle::close_status::exit_called == close_status && 0 == exit_status))
			{
				std::cerr << "ERROR: Subprocess with PID " << std::get <2>(res) << " exited with status " << exit_status;
				
				switch (close_status)
				{
					case lb::process_handle::close_status::unknown:
						std::cerr << " (exiting reason not known)";
						break;
					case lb::process_handle::close_status::terminated_by_signal:
						std::cerr << " (terminated by signal)";
						break;
					case lb::process_handle::close_status::stopped_by_signal:
						std::cerr << " (stopped by signal)";
						break;
					default:
						break;
				}
			
				std::cerr << '\n';
				std::exit(EXIT_FAILURE);
			}
		}
		else
		{
			lb::file_handle fh(lb::open_file_for_writing(dst_name, lb::writing_open_mode::CREATE));
			output_sequence_(ref_seq, graph, sample_idx, chr_copy_idx, fh);
		}
	}
	
	
	void output_sequences(sequence_vector const &ref_seq, variant_graph const &graph, char const * const pipe_cmd)
	{
		typedef variant_graph::ploidy_type ploidy_type;
		
		output_sequence(ref_seq, graph, variant_graph::SAMPLE_MAX, 0, pipe_cmd, "REF");
		for (auto const &[sample_idx, sample] : rsv::enumerate(graph.sample_names))
		{
			auto const ploidy(graph.sample_ploidy(sample_idx));
			for (auto const chr_copy_idx : rsv::iota(ploidy_type(0), ploidy))
			{
				// FIXME: Use std::format.
				std::stringstream dst_name;
				dst_name << sample;
				dst_name << '-';
				dst_name << chr_copy_idx;
				output_sequence(ref_seq, graph, sample_idx, chr_copy_idx, pipe_cmd, dst_name.str().data());
			}
		}
	}
	
	
	void run(char const *reference_path, char const *variants_path, char const *ref_seq_id, char const *chr_id, bool const should_output_sequences, char const *pipe_cmd, char const *graphviz_output_path)
	{
		// Read the reference sequence.
		sequence_vector ref_seq;
		{
			if (ref_seq_id)
				lb::log_time(std::cerr) << "Reading reference sequence with identifier “" << ref_seq_id << "”…" << std::flush;
			else
				lb::log_time(std::cerr) << "Reading the first reference sequence from the input FASTA…" << std::flush;
			auto const res(lb::read_single_fasta_sequence(reference_path, ref_seq, ref_seq_id));
			
			if (!res)
			{
				std::cerr << " ERROR: Unable to read the reference sequence.\n";
				std::exit(EXIT_FAILURE);
			}
			
			std::cerr << " Done. Reference length is " << ref_seq.size() << ".\n";
		}
		
		lb::log_time(std::cerr) << "Building the variant graph…" << std::flush;
		variant_graph graph;
		build_graph_statistics stats;
		build_variant_graph(ref_seq, variants_path, chr_id, graph, stats);
		std::cerr << " Done. Handled variants: " << stats.handled_variants << " chromosome ID mismatches: " << stats.chr_id_mismatches << "\n";
		
		if (graphviz_output_path)
		{
			lb::log_time(std::cerr) << "Outputting the variant graph in Grapnviz format…" << std::flush;
			lb::file_handle fh(lb::open_file_for_writing(graphviz_output_path, lb::writing_open_mode::CREATE));
			output_graphviz(ref_seq, graph, fh);
			std::cerr << " Done.\n";
		}
		
		if (should_output_sequences)
		{
			lb::log_time(std::cerr) << "Outputting sequences…" << std::flush;
			output_sequences(ref_seq, graph, pipe_cmd);
		}
	}
}


int main(int argc, char **argv)
{
#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif
	
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	if (args_info.show_invocation_given)
	{
		std::cerr << "Invocation:";
		for (int i(0); i < argc; ++i)
			std::cerr << ' ' << argv[i];
		std::cerr << '\n';
	}
	
	run(
		args_info.reference_arg,
		args_info.variants_arg,
		args_info.reference_sequence_arg,
		args_info.chromosome_arg,
		args_info.output_sequences_flag,
		args_info.pipe_arg,
		args_info.output_graphviz_arg
	);
	
	return EXIT_SUCCESS;
}
