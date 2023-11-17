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
	
	
	struct output_delegate
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
		output_delegate	*m_delegate{};
		
	public:
		output(char const *pipe_cmd, output_delegate &delegate):
			m_pipe_cmd(pipe_cmd),
			m_delegate(&delegate)
		{
		}
		
		virtual ~output() {}
		virtual void output_separate(sequence_type const &ref_seq, variant_graph const &graph) = 0;
		
		void output_a2m(sequence_type const &ref_seq, variant_graph const &graph, char const * const dst_name);
		
	protected:
		virtual void output_sequences_a2m(sequence_type const &ref_seq, variant_graph const &graph, libbio::file_ostream &stream) = 0;
		void output_sequence_file(sequence_type const &ref_seq, variant_graph const &graph, char const * const dst_name, sequence_writing_delegate &delegate);
	};
	
	
	class haplotype_output final : public output
	{
	public:
		using output::output;
		
		void output_separate(sequence_type const &ref_seq, variant_graph const &graph) override;
		
	protected:
		void output_sequences_a2m(sequence_type const &ref_seq, variant_graph const &graph, libbio::file_ostream &stream) override;
	};
	
	
	class founder_sequence_greedy_output final : public output
	{
	public:
		typedef variant_graph::ploidy_type					ploidy_type;
		typedef libbio::matrix <ploidy_type>				ploidy_matrix;
		typedef std::vector <variant_graph::position_type>	cut_position_vector;
		
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
		
	public:
		using output::output;
		
		void load_cut_positions(char const *path);
		void output_cut_positions(char const *path);
		[[nodiscard]] bool find_cut_positions(variant_graph const &graph, variant_graph::position_type const minimum_distance);
		
		[[nodiscard]] bool find_matching(variant_graph const &graph, ploidy_type const founder_count);
		
		void output_separate(sequence_type const &ref_seq, variant_graph const &graph) override;
		
	protected:
		void output_sequences_a2m(sequence_type const &ref_seq, variant_graph const &graph, libbio::file_ostream &stream) override;
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
