/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/adaptors.hpp>
#include <boost/range/combine.hpp>
#include <boost/container/map.hpp> // For an extract-capable map.
#include <vcf2multialign/tasks/read_subgraph_variants_task.hh>


namespace bm = boost::bimaps;


namespace vcf2multialign {
	
	void read_subgraph_variants_task::handle_variant(variant &var)
	{
		// Fill m_sequences_by_variant by pairs of sample_id -> sequence of alt indices.
		
		auto const var_pos(var.zero_based_pos());
		auto const &var_ref(var.ref());
		auto const var_ref_size(var_ref.size());
		auto const var_end(var_pos + var_ref_size);
		std::size_t i(0);
		m_variant_handler.enumerate_sample_genotypes(
			var,
			[
				this,
				&i,
				&var,
				var_end
			]
			(
				std::size_t const sample_no,
				uint8_t chr_idx,
				uint8_t const alt_idx,
				bool const is_phased
			) {
				always_assert(0 == chr_idx || is_phased, "Variant file not phased");
				always_assert(alt_idx <= 1); // FIXME: the maximum should be determined by checking the file before this step. Make use of int_vector from https://github.com/xxsds/sdsl-lite and use ceil(log2(max_chr_idx)) for the size.
				
				// If there is a variant, check that it does not overlap
				// with the co-ordinates of a previous variant in
				// the current sample. If it does, substitute with
				// a zero.
				if (chr_idx)
				{
					auto const var_pos(var.pos());
					if (var_pos < m_endpoints[i])
					{
						chr_idx = 0;

						// FIXME: log overlapping (per-sample) variants.
					}
					else
					{
						m_endpoints[i] = var_end;
					}
				}
				
				sample_id const sid(sample_no, chr_idx);
				auto it(m_sequences_by_sample.find(sid));
				if (m_sequences_by_sample.cend() == it)
				{
					auto &vec(m_sequences_by_sample[sid]);
					vec.reserve(m_variant_count);
					vec.emplace_back(alt_idx);
				}
				else
				{
					it->second.emplace_back(alt_idx);
				}
				
				++i;
			}
		);
	}
	
	
	void read_subgraph_variants_task::invert_sequences_by_sample(
		reduced_subgraph::sequence_vec &sequences,
		reduced_subgraph::sample_bimap &samples_by_sequence_idx
	)
	{
		// Use an std::map of sets instead of an std::multimap b.c. it makes moving the keys easier.
		typedef boost::container::map <
			reduced_subgraph::sequence_type,
			reduced_subgraph::sample_id_set
		> sample_map;
			
		typedef reduced_subgraph::sample_bimap::value_type sample_bimap_value;
		
		// Invert m_sequences_by_sample.
		// This step is needed mainly to keep only unique (ALT index) sequences.
		sample_map samples_by_sequence;
		while (m_sequences_by_sample.size())
		{
			auto it(m_sequences_by_sample.begin());
			auto &vec(it->second);
			samples_by_sequence[std::move(vec)].emplace(it->first);
		}
		
		std::size_t const unique_sequence_count(samples_by_sequence.size());
		
		// Number the vectors arbitrarily.
		// Reserve space for the sequences.
		sequences.clear();
		sequences.reserve(unique_sequence_count);
		
		// Move the sequences.
		while (samples_by_sequence.size())
		{
			auto it(samples_by_sequence.cbegin());
			auto node(samples_by_sequence.extract(it));
			
			auto &seq(node.key());
			auto &sample_ids(node.mapped());
			sequences.emplace_back(std::move(seq));
			
			auto const seq_idx(sequences.size() - 1);
			always_assert(seq_idx <= std::numeric_limits <reduced_subgraph::sequence_index>::max());
			
			for (auto const sample_id : sample_ids)
				samples_by_sequence_idx.insert(sample_bimap_value(seq_idx, sample_id));
		}
	}
	
	
	void read_subgraph_variants_task::split_sequences_to_paths(
		std::size_t const original_seq_count,
		reduced_subgraph::sample_bimap const &samples_by_sequence_idx,
		reduced_subgraph::path_map &generated_paths,
		reduced_subgraph::path_eq_map &generated_paths_eq
		//path_bimap &generated_paths,
		//path_eq_bimap &generated_paths_eq
	)
	{
		// Make sure that there are enough generated sequences.
		// FIXME: check this earlier to make error handling nicer? Or throw an exception if GCD allows that? Or call a delegate function?
		if (m_generated_path_count < original_seq_count)
		{
			std::cerr << original_seq_count << " generated sequences are needed, but the limit was " << m_generated_path_count << '.' << std::endl;
			abort();
		}
		
		auto const original_sample_count(samples_by_sequence_idx.right.size());
		std::vector <uint16_t> path_assignment_counts(original_seq_count, 0);	// In how many paths does a sequence of ALT indices appear.
		
		// Divide the generated paths by the ration of samples that have a particular path and the total sample count.
		{
			auto remaining_path_count(m_generated_path_count - original_seq_count);
			std::multimap <float, reduced_subgraph::sequence_index> fractionals;
			
			auto it(samples_by_sequence_idx.left.begin());
			auto const end(samples_by_sequence_idx.left.end());
			while (it != end)
			{
				// seq_id is an original sequence (vector of ALT indices) index.
				auto const seq_idx(it->first);
				always_assert(seq_idx < original_seq_count);
				
				// Unfortunately the iterator (likely) is not random access, so this will take linear time
				// (resulting in O(n) time complexity in total). lhs is decltype(seq_idx).
				auto const next(std::upper_bound(
					it, end, seq_idx, [](auto const lhs, auto const &rhs) -> bool { return lhs < rhs.first; }
				));
				
				// This could be done as part of finding next.
				auto const sample_count(std::distance(it, next));
				
				float const sample_count_fl(sample_count);
				float const ration(sample_count_fl / original_sample_count);
				
				std::size_t const integral(std::floor(ration));
				float const fractional(ration - integral);
				
				assert(integral <= remaining_path_count);
			
				// Associate the fractional part with the sequence index.
				fractionals.emplace(
					std::piecewise_construct,
					std::forward_as_tuple(fractional),
					std::forward_as_tuple(seq_idx)
				);
			
				// Store the integral part. Make sure that at least one sequence is assigned.
				auto const count(1 + integral);
				{
					typedef decltype(path_assignment_counts)::value_type assignment_count_type;
					always_assert(count < std::numeric_limits <assignment_count_type>::max());
				}
				path_assignment_counts[seq_idx] = count;
			
				// Update remaining_path_count;
				remaining_path_count -= integral;
				
				it = std::move(next);
			}
			
			// Add generated sequences in the order of fractionals.
			for (auto const &kv : boost::adaptors::reverse(fractionals))
			{
				// kv.first is the fractional part.
				// kv.second is an original sequence (vector of ALT indices) index.
				auto const seq_idx(kv.second);
				if (0 == remaining_path_count)
					goto loop_end;
				
				++path_assignment_counts[seq_idx];
				--remaining_path_count;
			}
			
		loop_end:
			always_assert(0 == remaining_path_count);
		}
		
		// Split the original sequence identifiers among paths by the assignment counts.
		{
			std::size_t seq_idx(0);
			std::size_t path_idx(0);
			for (auto const &count : path_assignment_counts)
			{
				always_assert(seq_idx <= std::numeric_limits <reduced_subgraph::sequence_index>::max());
				auto current_remaining_count(count);
				auto const range(samples_by_sequence_idx.left.equal_range(seq_idx)); // FIXME: this is O(log n) as long as samples_by_sequence_idx.left is a set. If changing the type does not help, find the next endpoint of the range with std::upper_bound to get O(n) time.
				auto it(range.first);
				auto const end(range.second);
				
				// Store the equivalence classes first.
				for (std::size_t i(0); i < count; ++i)
				{
					auto const current_path_idx(i + path_idx);
					always_assert(current_path_idx <= std::numeric_limits <reduced_subgraph::path_index>::max());
					//generated_paths_eq.insert(current_path_idx, seq_idx);
					generated_paths_eq.emplace(
						std::piecewise_construct,
						std::forward_as_tuple(current_path_idx),
						std::forward_as_tuple(seq_idx)
					);
				}
				
				// Map the sequence identifiers to paths.
				bool check(false);
				while (it != end)
				{
					for (std::size_t i(0); i < count; ++i)
					{
						//generated_paths.insert(i + path_idx, *it);
						generated_paths.emplace(
							std::piecewise_construct,
							std::forward_as_tuple(i + path_idx),
							std::forward_as_tuple(it->second)
						);
						++it;
						
						if (end == it)
							goto loop_end_2;
					}
					
					// Check that all the paths have a sequence identifier.
					check = true;
				}
				
			loop_end_2:
				always_assert(check);
				path_idx += count;
				++seq_idx;
			}
		}
	}
	
	
	void read_subgraph_variants_task::finish()
	{
		reduced_subgraph::sequence_vec	sequences;
		reduced_subgraph::sample_bimap	samples_by_sequence_idx;
		reduced_subgraph::path_map		generated_paths;
		reduced_subgraph::path_eq_map	generated_paths_eq;
		
		// Finished processing the variants.
		// Use the variant sequences as keys.
		invert_sequences_by_sample(sequences, samples_by_sequence_idx);
		// m_samples_by_sequence now contains variant sequences as keys and sample identifiers as values.
		
		split_sequences_to_paths(
			sequences.size(),
			samples_by_sequence_idx,
			generated_paths,
			generated_paths_eq
		);
			
		reduced_subgraph rsg(
			std::move(sequences),
			std::move(samples_by_sequence_idx),
			std::move(generated_paths),
			std::move(generated_paths_eq),
			m_start_lineno,
			m_variant_count
		);
		
		m_delegate->task_did_finish(*this, std::move(rsg));
	}
	
	
	void read_subgraph_variants_task::execute()
	{
		m_vcf_reader.set_parsed_fields(vcf_field::ALL);
		m_variant_handler.process_variants();
	}
}
