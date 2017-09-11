/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/container/list.hpp>
#include <list>
#include <vcf2multialign/sample_reducer.hh>


namespace vcf2multialign {
	
	std::ostream &operator<<(std::ostream &stream, variant_sequence const &seq)
	{
		stream
		<< "Start: " << seq.start_pos()
		<< " end: " << seq.end_pos()
		<< " length: " << seq.length()
		<< " sample: " << seq.sample_no()
		<< " chr: " << (int) seq.chr_idx()
		<< " sequence:";

		for (auto const &kv : seq.m_alt_indices)
			stream << " (" << kv.first << ", " << (int) kv.second << ")";

		return stream;
	}
	
	
	// Check whether prepared_sequences already contains seq.
	bool sample_reducer::prepared_contains_sequence(variant_sequence const &seq) const
	{
		auto it(m_prepared_sequences.find(seq.start_pos()));
		if (m_prepared_sequences.cend() == it)
			return false;
		
		for (auto const &other_seq : it->second)
		{
			if (seq.equal_sequences(other_seq))
				return true;
		}
		
		return false;
	}
	
	
	// Move seq to prepared_sequences if it hasn't been already added.
	void sample_reducer::check_and_copy_seq_to_prepared(variant_sequence &seq)
	{
		auto const start_pos(seq.start_pos_1());
		if (!(0 == start_pos || prepared_contains_sequence(seq)))
		{
			auto &vec(m_prepared_sequences[start_pos - 1]);
			auto &dst(vec.emplace_back(seq.seq_id()));
			
			using std::swap;
			swap(dst, seq);
		}
	}
	
	
	// Move variant_seq to prepared_sequences if there are no variants in the padding distance.
	bool sample_reducer::check_variant_sequence(
		variant_sequence &seq,
		variant_sequence_id const &seq_id,
		std::size_t const zero_based_pos
	)
	{
		if (seq.assign_id(seq_id))
		{
			seq.set_start_pos(zero_based_pos);
			return true;
		}
		
		auto const end_pos(seq.end_pos());
		if (! (end_pos <= zero_based_pos))
			return false;
		
		if (m_padding_amt < zero_based_pos - end_pos)
		{
			check_and_copy_seq_to_prepared(seq);
			seq.reset();
			seq.set_start_pos(zero_based_pos);
		}
		
		return true;
	}
	
	
	// Create subsequences.
	void sample_reducer::handle_variant(variant &var)
	{
		// Verify that the positions are in increasing order.
		auto const lineno(var.lineno());
		auto const pos(var.zero_based_pos());
		auto const end(pos + var.ref().size());
		
		always_assert(m_last_position <= pos, "Positions not in increasing order");
		
		for (auto const &kv : m_delegate->sample_names())
		{
			auto const sample_no(kv.second);
			
			// Handle the genotype.
			m_delegate->enumerate_genotype(var, sample_no,
				[
					this,
					lineno,
					pos,
					end,
					sample_no
				](
					uint8_t const chr_idx, std::size_t const alt_idx, bool const is_phased
				) {
					always_assert(0 == chr_idx || is_phased, "Variant file not phased");
			
					if (0 != alt_idx && m_delegate->is_valid_alt(alt_idx))
					{
						// If the current variant references an ALT column value in a chromosome, add the
						// ALT index to the corresponding sequence.
						if (alt_idx)
						{
							variant_sequence_id seq_id(sample_no, chr_idx);
							variant_sequence &seq(m_variant_sequences[seq_id]);
					
							// First check if the previous variant is beyond the padding distance.
							if (check_variant_sequence(seq, seq_id, pos))
							{
								seq.add_alt(lineno, pos, alt_idx);
								m_delegate->assigned_alt_to_sequence(alt_idx);
							}
							else
							{
								auto const sample_no(seq.sample_no());
								auto const chr_idx(seq.chr_idx());
								m_delegate->found_overlapping_alt(lineno, alt_idx, sample_no, chr_idx);
							}
					
							m_delegate->handled_alt(alt_idx);
						}
					}
				}
			);
		}
			
		m_last_position = pos;
	}
	
	
	void sample_reducer::prepare()
	{
		m_last_position = 0;
	}
	
	
	void sample_reducer::finish()
	{
		// Add the remaining ranges.
		for (auto &kv : m_variant_sequences)
		{
			auto &seq(kv.second);
			check_and_copy_seq_to_prepared(seq);
		}
		
		assign_ranges_greedy();
	}
	
	
	void sample_reducer::print_prepared_sequences()
	{
		for (auto const &kv : m_prepared_sequences)
		{
			auto const key(kv.first);
			auto const &list(kv.second);
			std::cerr << "Key: " << key << std::endl;
			for (auto const &item : list)
				std::cerr << "\t" << item << std::endl;
		}
	}
	
	
	void sample_reducer::assign_ranges_greedy()
	{
		//print_prepared_sequences(prepared_sequences);

		// Sort each range vector.
		for (auto &kv : m_prepared_sequences)
		{
			auto &vec(kv.second);
			vec.sort([](variant_sequence const &a, variant_sequence const &b) {
				return (a.length() > b.length());
			});
		}
		
		//std::cerr << "Sorted the range vectors." << std::endl;
		//print_prepared_sequences(prepared_sequences);

		std::size_t compressed_idx(0);

		// If reference is to be output, add an empty compressed range for it.
		static_assert(0 == REF_SAMPLE_NUMBER);
		if (m_output_ref)
		{
			m_compressed_ranges->emplace_back();
			++compressed_idx;
		}
		
		// Assign ranges to sequence identifiers in a greedy manner.
		// prepared_sequences is indexed with 1-based start positions while
		while (m_prepared_sequences.size())
		{
			std::size_t end_pos(0);
			
			// For some reason retrieving reference from emplace_back causes problems.
			m_compressed_ranges->emplace_back();
			auto &dst(*m_compressed_ranges->rbegin());
			
			//std::cerr << "Assigning the following sequences to compressed sequence " << compressed_idx << ':' << std::endl;
			while (true)
			{
				// Find the list which contains ranges that start from end_pos.
				auto next_it(m_prepared_sequences.lower_bound(end_pos));
				if (m_prepared_sequences.cend() == next_it)
					break;
				
				auto &list(next_it->second);
				assert(list.size());
				auto const next_seq_it(list.cbegin());
				auto const start_pos(next_seq_it->start_pos());
				end_pos = next_seq_it->end_pos() + m_padding_amt;
				
				// Move the range to the current list.
				//std::cerr << '\t' << (*next_seq_it) << std::endl;
				dst.emplace_hint(dst.cend(), start_pos, std::move(*next_seq_it));
				
				list.erase(next_seq_it);
				if (0 == list.size())
					m_prepared_sequences.erase(next_it);
			}
			
			++compressed_idx;
		}
	}
	
	
#if 0
	void reduce_samples(
		vcf_reader &reader,
		error_logger &error_logger,
		variant_set const &skipped_variants,
		v2m::sv_handling const sv_handling,
		std::size_t const padding_amt,
		bool const output_ref,
		range_map &compressed_ranges
	)
	{
		sample_reducer compressor(sv_handling, error_logger);
		
		subsequence_map prepared_sequences;
		
		std::cerr << "Compressing variants… " << std::endl;
		create_subsequences(reader, error_logger, skipped_variants, padding_amt, prepared_sequences);
		
		std::cerr << "Assigning variant ranges to new haplotype sequences… " << std::flush;
		assign_ranges_greedy(prepared_sequences, compressed_ranges, padding_amt, output_ref);
		std::cerr << "expressing variants with " << compressed_ranges.size() << " sequences." << std::endl;
	}
#endif
}
