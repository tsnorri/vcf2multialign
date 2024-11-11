/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_HH

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <libbio/assert.hh>
#include <libbio/file_handling.hh>
#include <libbio/int_matrix.hh>
#include <libbio/int_matrix/cereal_serialization.hh>
#include <libbio/int_matrix/size_calculation.hh>
#include <libbio/int_vector/cereal_serialization.hh>
#include <libbio/int_vector/size_calculation.hh>
#include <libbio/size_calculator.hh>
#include <libbio/vcf/variant.hh>
#include <limits>										// std::numeric_limits
#include <range/v3/view/zip.hpp>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>										// std::pair
#include <vector>


namespace vcf2multialign {

	typedef std::uint32_t					cereal_version_type;
	typedef std::vector <char>				sequence_type;


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
		typedef libbio::bit_matrix			path_matrix;

		constexpr static inline position_type const POSITION_MAX{std::numeric_limits <position_type>::max()};
		constexpr static inline node_type const NODE_MAX{std::numeric_limits <node_type>::max()};
		constexpr static inline edge_type const EDGE_MAX{std::numeric_limits <edge_type>::max()};
		constexpr static inline sample_type const SAMPLE_MAX{std::numeric_limits <sample_type>::max()};
		constexpr static inline ploidy_type const PLOIDY_MAX{std::numeric_limits <ploidy_type>::max()};

		position_vector						reference_positions;			// Reference positions by node number.
		position_vector						aligned_positions;				// MSA positions by node number.
		node_vector							alt_edge_targets;				// ALT edge targets by edge number.
		edge_vector							alt_edge_count_csum;			// Cumulative sum of ALT edge counts by 1-based node number.
		label_vector						alt_edge_labels;				// ALT edge labels by edge number.
		path_matrix							paths_by_chrom_copy_and_edge;	// Edges on rows, chromosome copies (samples multiplied by ploidy) in columns.
		path_matrix							paths_by_edge_and_chrom_copy;	// Chromosome copies on rows, edges in columns.

		label_vector						sample_names;					// Sample names by sample index. FIXME: In case we have variant_graph ->> chromosome at some point, this should be in the graph.
		ploidy_csum_vector					ploidy_csum;					// Cumulative sum of ploidies by 1-based sample number (for this chromosome).

		node_type node_count() const { return reference_positions.size(); }
		edge_type edge_count() const { return alt_edge_targets.size(); }

		std::pair <edge_type, edge_type> edge_range_for_node(node_type const &node_idx) const { return {alt_edge_count_csum[node_idx], alt_edge_count_csum[1 + node_idx]}; }

		ploidy_type sample_ploidy(sample_type const sample_idx) const { return ploidy_csum[1 + sample_idx] - ploidy_csum[sample_idx]; }
		ploidy_type total_chromosome_copies() const { return ploidy_csum.back(); }

		node_type add_node(position_type const ref_pos, position_type const aln_pos);
		node_type add_or_update_node(position_type const ref_pos, position_type const aln_pos);
		edge_type add_edge(std::string_view const label = std::string_view{});

		position_type aligned_length(node_type const lhs, node_type const rhs) const { return aligned_positions[rhs] - aligned_positions[lhs]; }

		// For Cereal
		template <typename t_archive> void serialize(t_archive &ar, cereal_version_type const version);
	};


	class variant_graph_walker
	{
	public:
		typedef variant_graph::node_type		node_type;
		typedef variant_graph::edge_type		edge_type;
		typedef variant_graph::position_type	position_type;

		static_assert(std::is_unsigned_v <node_type>);

	protected:
		sequence_type const	*m_reference{};
		variant_graph const	*m_graph{};
		node_type			m_node{variant_graph::NODE_MAX};

	public:
		variant_graph_walker() = default;

		explicit variant_graph_walker(variant_graph const &graph):
			m_graph(&graph)
		{
		}

		variant_graph_walker(sequence_type const &reference, variant_graph const &graph):
			m_reference(&reference),
			m_graph(&graph)
		{
		}

		bool advance() { return ++m_node < m_graph->node_count(); }
		variant_graph const &graph() const { return *m_graph; }
		node_type node() const { return m_node; }
		position_type ref_position() const { return m_graph->reference_positions[m_node]; }
		position_type aligned_position() const { return m_graph->aligned_positions[m_node]; }
		std::string_view ref_label() const { return ref_label_(m_node); }
		std::string_view ref_label(std::size_t const rhs_node) const { libbio_assert_lte(m_node, rhs_node); return ref_label_(rhs_node); }
		position_type ref_length() const { return ref_length_(m_node); }
		position_type ref_length(std::size_t const rhs_node) const { libbio_assert_lte(m_node, rhs_node); return ref_length_(rhs_node); }
		position_type aligned_length() const { return aligned_length_(m_node); }
		position_type aligned_length(std::size_t const rhs_node) const { libbio_assert_lte(m_node, rhs_node); return aligned_length_(rhs_node); }
		edge_type alt_edge_count() const { return m_graph->alt_edge_count_csum[1 + m_node] - m_graph->alt_edge_count_csum[m_node]; }
		edge_type alt_edge_base() const { return m_graph->alt_edge_count_csum[m_node]; }
		inline auto alt_edge_labels() const;
		inline auto alt_edge_targets() const;
		inline auto alt_edges() const;

	protected:
		std::string_view ref_label_(std::size_t const rhs_node) const { auto const *data(m_reference->data()); return {data + ref_position(), data + m_graph->reference_positions[rhs_node]}; }
		inline std::size_t ref_length_(std::size_t const rhs_node) const { return m_graph->reference_positions[rhs_node] - ref_position(); }
		inline std::size_t aligned_length_(std::size_t const rhs_node) const { return m_graph->aligned_positions[rhs_node] - aligned_position(); }
	};


	struct build_graph_delegate
	{
		typedef variant_graph::position_type	position_type;

		virtual ~build_graph_delegate() {}
		virtual bool should_include(std::string_view const sample_name, variant_graph::ploidy_type const chrom_copy_idx) const = 0;

		// FIXME: add typedef for line number.
		virtual void report_overlapping_alternative(
			std::uint64_t const lineno,
			variant_graph::position_type const ref_pos,
			std::vector <std::string_view> const &var_id,
			std::string_view const sample_name,
			variant_graph::ploidy_type const chrom_copy_idx,
			std::uint32_t const gt
		) = 0;

		// FIXME: add typedef for variant index.
		virtual bool ref_column_mismatch(std::uint64_t const var_idx, libbio::vcf::transient_variant const &var, std::string_view const expected) = 0;
	};


	struct process_graph_delegate
	{
		virtual ~process_graph_delegate() {}
		virtual void handled_node(variant_graph::node_type const node) = 0;
	};


	struct build_graph_statistics
	{
		std::uint64_t	handled_variants{};
		std::uint64_t	chr_id_mismatches{};
	};


	void build_variant_graph(
		sequence_type const &ref_seq,
		char const *variants_path,
		char const *chr_id,
		variant_graph &graph,
		build_graph_statistics &stats,
		build_graph_delegate &delegate
	);


	inline void build_variant_graph(
		sequence_type const &ref_seq,
		std::filesystem::path const &variants_path,
		char const *chr_id,
		variant_graph &graph,
		build_graph_statistics &stats,
		build_graph_delegate &delegate
	)
	{
		build_variant_graph(ref_seq, variants_path.c_str(), chr_id, graph, stats, delegate);
	}


	void phase(variant_graph &graph, std::uint16_t const ploidy, libbio::file_ostream &os);


	template <typename t_archive>
	void variant_graph::serialize(t_archive &ar, cereal_version_type const version)
	{
		ar(reference_positions);
		ar(aligned_positions);
		ar(alt_edge_targets);
		ar(alt_edge_count_csum);
		ar(alt_edge_labels);
		ar(paths_by_chrom_copy_and_edge);
		ar(paths_by_edge_and_chrom_copy);
		ar(sample_names);
		ar(ploidy_csum);
	}


	auto variant_graph_walker::alt_edge_labels() const
	{
		auto const *data(m_graph->alt_edge_labels.data());
		return std::span(data + m_graph->alt_edge_count_csum[m_node], data + m_graph->alt_edge_count_csum[1 + m_node]);
	}


	auto variant_graph_walker::alt_edge_targets() const
	{
		auto const *data(m_graph->alt_edge_targets.data());
		return std::span(data + m_graph->alt_edge_count_csum[m_node], data + m_graph->alt_edge_count_csum[1 + m_node]);
	}


	auto variant_graph_walker::alt_edges() const
	{
		return ranges::views::zip(alt_edge_targets(), alt_edge_labels());
	}
}


namespace libbio::size_calculation {

	template <>
	struct value_size_calculator <vcf2multialign::variant_graph>
	{
		void operator()(size_calculator &sc, entry_index_type const entry_idx, vcf2multialign::variant_graph const &graph) const;
	};
}

#endif
