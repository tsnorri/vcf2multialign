/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_OUTPUT_HH
#define VCF2MULTIALIGN_OUTPUT_HH

#include <libbio/matrix.hh>
#include <libbio/subprocess.hh>
#include <vcf2multialign/find_cut_positions.hh>
#include <vcf2multialign/variant_graph.hh>


namespace vcf2multialign {
	
	typedef libbio::subprocess <libbio::subprocess_handle_spec::STDIN>	subprocess_type;
	
	
	struct sequence_writing_delegate; // Fwd.
	
	
	struct output_delegate : public process_graph_delegate
	{
		typedef variant_graph::sample_type	sample_type;
		typedef variant_graph::ploidy_type	ploidy_type;
		typedef std::uint32_t				sequence_count_type;
		
		virtual ~output_delegate() {}
		virtual void will_handle_sample(std::string const &sample, sample_type const sample_idx, ploidy_type const chr_copy_idx) = 0;
		virtual void will_handle_founder_sequence(sample_type const idx) = 0;
		virtual void handled_sequences(sequence_count_type const sequence_count) = 0;
		virtual void exit_subprocess(subprocess_type &proc) = 0;
	};
	
	
	class output
	{
	protected:
		char const		*m_pipe_cmd{};
		char const		*m_chromosome_id{};
		output_delegate	*m_delegate{};
		bool			m_should_output_reference{};
		bool			m_should_output_unaligned{};
		
	public:
		output(char const *pipe_cmd, char const *chromosome_id, bool const should_output_reference, bool const should_output_unaligned, output_delegate &delegate):
			m_pipe_cmd(pipe_cmd),
			m_chromosome_id(chromosome_id),
			m_delegate(&delegate),
			m_should_output_reference(should_output_reference),
			m_should_output_unaligned(should_output_unaligned)
		{
		}
		
		virtual ~output() {}
		virtual void output_separate(sequence_type const &ref_seq, variant_graph const &graph, bool const should_include_fasta_header) = 0;
		
		void output_a2m(sequence_type const &ref_seq, variant_graph const &graph, char const * const dst_name);
		virtual void output_a2m(sequence_type const &ref_seq, variant_graph const &graph, std::ostream &stream) = 0;
		
	protected:
		void output_sequence_file(sequence_type const &ref_seq, variant_graph const &graph, char const * const dst_name, bool const should_include_fasta_header, sequence_writing_delegate &delegate);
	};
	
	
	class haplotype_output final : public output
	{
	public:
		using output::output;
		
		void output_separate(sequence_type const &ref_seq, variant_graph const &graph, bool const should_include_fasta_header) override;
		void output_a2m(sequence_type const &ref_seq, variant_graph const &graph, std::ostream &stream) override;
	};
	
	
	class founder_sequence_greedy_output final : public output
	{
	public:
		typedef variant_graph::ploidy_type					ploidy_type;
		typedef libbio::matrix <ploidy_type>				ploidy_matrix;
		
		constexpr static inline auto const PLOIDY_MAX{variant_graph::PLOIDY_MAX};
		
		struct cut_positions
		{
			cut_position_vector			cut_positions;
			variant_graph::edge_type	min_distance{};
			cut_position_score_type		score{};
			
			// For Cereal.
			template <typename t_archive> void serialize(t_archive &ar, cereal_version_type const version);
		};
		
	private:
		cut_positions	m_cut_positions;
		ploidy_matrix	m_assigned_samples;
		bool			m_should_keep_ref_edges{};
		
	public:
		founder_sequence_greedy_output(
			char const *pipe_cmd,
			char const *chromosome_id,
			bool const should_output_reference,
			bool const should_keep_ref_edges,
			bool const should_output_unaligned,
			output_delegate &delegate
		):
			output(pipe_cmd, chromosome_id, should_output_reference, should_output_unaligned, delegate),
			m_should_keep_ref_edges(should_keep_ref_edges)
		{
		}
		
		[[nodiscard]] cut_position_vector const &cut_positions() const { return m_cut_positions.cut_positions; }
		[[nodiscard]] ploidy_matrix const &assigned_samples() const { return m_assigned_samples; }

		void load_cut_positions(char const *path);
		void output_cut_positions(char const *path);
		[[nodiscard]] bool find_cut_positions(variant_graph const &graph, variant_graph::position_type const minimum_distance);
		[[nodiscard]] cut_position_score_type max_segmentation_height() const { return m_cut_positions.score; }
		
		[[nodiscard]] bool find_matchings(variant_graph const &graph, ploidy_type const founder_count);
		
		void output_separate(sequence_type const &ref_seq, variant_graph const &graph, bool const should_include_fasta_header) override;
		void output_a2m(sequence_type const &ref_seq, variant_graph const &graph, std::ostream &stream) override;
	};
	
	
	template <typename t_archive>
	void founder_sequence_greedy_output::cut_positions::serialize(t_archive &ar, cereal_version_type const version)
	{
		ar(min_distance);
		ar(cut_positions);
		ar(score);
	}
}

#endif
