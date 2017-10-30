/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_REDUCED_SUBGRAPH_HH
#define VCF2MULTIALIGN_REDUCED_SUBGRAPH_HH

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/format.hpp>
#include <map>
#include <sdsl/int_vector.hpp>
#include <set>
#include <vcf2multialign/util.hh>


namespace vcf2multialign {
	// Use 32 bits for sample ids.
	struct sample_id {
		enum
		{
			SAMPLE_BIT = 24,
			CHR_BIT = 8
		};
		
		uint32_t sample : SAMPLE_BIT;
		uint8_t chr : CHR_BIT;
		
		sample_id(std::size_t sample_, std::size_t chr_):
			sample(sample_),
			chr(chr_)
		{
			always_assert(sample_ < 1 << (SAMPLE_BIT - 1));
			always_assert(chr_ < 1 << (CHR_BIT - 1));
		}
		
		sample_id(): sample_id(0, 0) {}
		
		bool operator<(sample_id const &other) const
		{
			return std::tie(sample, chr) < std::tie(other.sample, other.chr);
		}
	};
	
	static_assert(4 == sizeof(sample_id));
	
	
	// FIXME: everything that uses consecutive indices here should be replaced with std::vectors or bimaps with vector views for O(1) random access.
	class reduced_subgraph
	{
	public:
		typedef sdsl::int_vector <0> sequence_type;
		
		typedef std::set <sample_id> sample_id_set;
		
		typedef uint16_t sequence_index;	// FIXME: could be uint8_t?
		typedef uint16_t path_index;		// FIXME: could be uint8_t?
	
		// Map unique sequence indices or path indices -> (ALT index) sequences.
		typedef std::vector <sequence_type> sequence_vec;
		
		// Map (existing) unique (ALT index) sequence indices <->> (existing) sample indices.
		// We use a bimap in order to avoid creating equivalence classes of sample ids (by paths).
		// Instead, we store the equivalence class (sequence index) of each path and use this container
		// to compare.
		typedef boost::bimap <
			boost::bimaps::set_of <sequence_index>,
			boost::bimaps::multiset_of <sample_id>
		> sample_bimap;
		
		// Map generated path indices <->> sample id.
		typedef std::multimap <path_index, sample_id> path_map;
	
		// Map generated path indices <<-> (existing) unique (ALT index) sequence index.
		typedef std::map <path_index, sequence_index> path_eq_map;
		
	protected:
		sequence_vec						m_sequences;				// Initialized in invert_sequences_by_sample().
		sample_bimap						m_samples_by_sequence_idx;	// Initialized in invert_sequences_by_sample().
		path_map							m_generated_paths;
		path_eq_map							m_generated_paths_eq;
		std::size_t							m_start_lineno{0};
		std::size_t							m_variant_count{0};
		
	public:
		reduced_subgraph() = default;
		
		reduced_subgraph(
			sequence_vec &&sequences,
			sample_bimap &&samples_by_sequence_idx,
			path_map &&generated_paths,
			path_eq_map &&generated_paths_eq,
			std::size_t const start_lineno,
			std::size_t const variant_count
		):
			m_sequences(std::move(sequences)),
			m_samples_by_sequence_idx(std::move(samples_by_sequence_idx)),
			m_generated_paths(std::move(generated_paths)),
			m_generated_paths_eq(std::move(generated_paths_eq)),
			m_start_lineno(start_lineno),
			m_variant_count(variant_count)
		{
		}
		
		std::size_t start_lineno() const { return m_start_lineno; }
		std::size_t variant_count() const { return m_variant_count; }
		
		sequence_type const &path_sequence(path_index const idx) const
		{
			sequence_index seq_idx{0};
			bool const st(path_sequence_index(idx, seq_idx));
			always_assert(st);
			return m_sequences.at(seq_idx);
		}
		
		std::pair <path_map::const_iterator, path_map::const_iterator>
		path_samples(path_index const idx) const
		{
			return m_generated_paths.equal_range(idx);
		}
		
		bool sample_sequence_index(sample_id const &sample_id, sequence_index /* out */ &seq_idx) const
		{
			auto const seq_it(m_samples_by_sequence_idx.right.find(sample_id));
			if (seq_it == m_samples_by_sequence_idx.right.end())
				return false;
			
			seq_idx = seq_it->second;
			return true;
		}
		
		bool path_sequence_index(path_index const &path_idx, sequence_index /* out */ &seq_idx) const
		{
			auto const seq_it(m_generated_paths_eq.find(path_idx));
			if (m_generated_paths_eq.cend() != seq_it)
				return false;
			
			seq_idx = seq_it->second;
			return true;
		}
	};
	
	std::size_t edge_weight(
		reduced_subgraph const &lhs,
		reduced_subgraph const &rhs,
		reduced_subgraph::path_index const li,
		reduced_subgraph::path_index const ri
	);
}

#endif
