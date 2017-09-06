/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/container/list.hpp>
#include <list>
#include <vcf2multialign/compress_variants.hh>
#include <vcf2multialign/variant_buffer.hh>


namespace vcf2multialign {
	
	typedef std::map <std::size_t, boost::container::list <variant_sequence>> subsequence_map;

	
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
	bool contains_sequence(
		subsequence_map const &prepared_sequences,
		variant_sequence const &seq
	)
	{
		auto it(prepared_sequences.find(seq.start_pos()));
		if (prepared_sequences.cend() == it)
			return false;
		
		for (auto const &other_seq : it->second)
		{
			if (seq.equal_sequences(other_seq))
				return true;
		}
		
		return false;
	}
	
	
	// Move seq to prepared_sequences if it hasn't been already added.
	void check_and_copy_seq_to_set(
		variant_sequence &seq,
		subsequence_map &prepared_sequences // std::map <std::size_t, boost::container::list <variant_sequence>> subsequence_map
	)
	{
		auto const start_pos(seq.start_pos_1());
		if (!(0 == start_pos || contains_sequence(prepared_sequences, seq)))
		{
			auto &vec(prepared_sequences[start_pos - 1]);
			auto &dst(vec.emplace_back(seq.seq_id()));
			
			using std::swap;
			swap(dst, seq);
		}
	}
	
	
	// Move variant_seq to prepared_sequences if there are no variants in the padding distance.
	void check_variant_sequence(
		variant_sequence &seq,
		variant_sequence_id const &seq_id,
		std::size_t const zero_based_pos,
		std::size_t const padding_amt,
		subsequence_map &prepared_sequences
	)
	{
		if (seq.assign_id(seq_id))
		{
			seq.set_start_pos(zero_based_pos);
			return;
		}
		
		auto const end_pos(seq.end_pos());
		assert(end_pos <= zero_based_pos);
		if (padding_amt < zero_based_pos - end_pos)
		{
			check_and_copy_seq_to_set(seq, prepared_sequences);
			seq.reset();
			seq.set_start_pos(zero_based_pos);
		}
	}
	
	
	void create_subsequences(
		vcf_reader &reader,
		error_logger &error_logger,
		variant_set const &skipped_variants,
		std::size_t const padding_amt,
		subsequence_map &prepared_sequences
	)
	{
		size_t last_position(0);
		size_t i(0);
		
		// Variant sequences by sample number.
		std::map <variant_sequence_id, variant_sequence> variant_sequences;
		
		reader.reset();
		reader.set_parsed_fields(vcf_field::ALL);
		bool should_continue(false);
		do {
			reader.fill_buffer();
			should_continue = reader.parse(
				[
					&reader,
					&error_logger,
					&skipped_variants,
					&padding_amt,
					&prepared_sequences,
					&last_position,
					&i,
					&variant_sequences
				]
				(transient_variant const &var)
				-> bool
				{
					auto const lineno(var.lineno());
					if (0 != skipped_variants.count(lineno))
						return true;
					
					// Verify that the positions are in increasing order.
					auto const pos(var.zero_based_pos());
					
					always_assert(last_position <= pos, "Positions not in increasing order");
					
					auto const end(pos + var.ref().size());
					auto const var_lineno(var.lineno());
					
					for (auto const &kv : reader.sample_names())
					{
						auto const sample_no(kv.second);
						auto const &sample(var.sample(sample_no));
					
						// Handle the genotype.
						uint8_t chr_idx(0);
						for (auto const gt : sample.get_genotype())
						{
							auto const alt_idx(gt.alt);
							auto const is_phased(gt.is_phased);
							always_assert(0 == chr_idx || is_phased, "Variant file not phased");
							
							// If the current variant references an ALT column value in a chromosome, add the
							// ALT index to the corresponding sequence.
							if (alt_idx)
							{
								variant_sequence_id seq_id(sample_no, chr_idx);
								variant_sequence &seq(variant_sequences[seq_id]);
								
								// First check if the previous variant is beyond the padding distance.
								check_variant_sequence(seq, seq_id, pos, padding_amt, prepared_sequences);
								
								seq.add_alt(var_lineno, pos, alt_idx);
							}
							
							++chr_idx;
						}
					}
						
					last_position = pos;
					
					++i;
					if (0 == i % 100000)
						std::cerr << "Handled " << i << " variants…" << std::endl;
				
					return true;
				}
			);
		} while (should_continue);
		
		// Add the remaining ranges.
		for (auto &kv : variant_sequences)
		{
			auto &seq(kv.second);
			check_and_copy_seq_to_set(seq, prepared_sequences);
		}
	}
	
	
	void print_prepared_sequences(subsequence_map const &prepared_sequences)
	{
		for (auto const &kv : prepared_sequences)
		{
			auto const key(kv.first);
			auto const &list(kv.second);
			std::cerr << "Key: " << key << std::endl;
			for (auto const &item : list)
				std::cerr << "\t" << item << std::endl;
		}
	}
	
	
	void assign_ranges_greedy(
		subsequence_map &prepared_sequences,
		range_map &compressed_ranges,
		std::size_t const padding_amt,
		bool const output_ref
	)
	{
		//print_prepared_sequences(prepared_sequences);

		// Sort each range vector.
		for (auto &kv : prepared_sequences)
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
		if (output_ref)
		{
			compressed_ranges.emplace_back();
			++compressed_idx;
		}
		
		// Assign ranges to sequence identifiers in a greedy manner.
		// prepared_sequences is indexed with 1-based start positions while
		while (prepared_sequences.size())
		{
			std::size_t end_pos(0);
			
			// For some reason retrieving reference from emplace_back causes problems.
			compressed_ranges.emplace_back();
			auto &dst(*compressed_ranges.rbegin());
			
			//std::cerr << "Assigning the following sequences to compressed sequence " << compressed_idx << ':' << std::endl;
			while (true)
			{
				// Find the list which contains ranges that start from end_pos.
				auto next_it(prepared_sequences.lower_bound(end_pos));
				if (prepared_sequences.cend() == next_it)
					break;
				
				auto &list(next_it->second);
				assert(list.size());
				auto const next_seq_it(list.cbegin());
				auto const start_pos(next_seq_it->start_pos());
				end_pos = next_seq_it->end_pos() + padding_amt;
				
				// Move the range to the current list.
				//std::cerr << '\t' << (*next_seq_it) << std::endl;
				dst.emplace_hint(dst.cend(), start_pos, std::move(*next_seq_it));
				
				list.erase(next_seq_it);
				if (0 == list.size())
					prepared_sequences.erase(next_it);
			}
			
			++compressed_idx;
		}
	}
	
	
	void compress_variants(
		vcf_reader &reader,
		error_logger &error_logger,
		variant_set const &skipped_variants,
		std::size_t const padding_amt,
		bool const output_ref,
		range_map &compressed_ranges
	)
	{
		subsequence_map prepared_sequences;
		
		std::cerr << "Compressing variants… " << std::endl;
		create_subsequences(reader, error_logger, skipped_variants, padding_amt, prepared_sequences);
		
		std::cerr << "Assigning variant ranges to new haplotype sequences… " << std::flush;
		assign_ranges_greedy(prepared_sequences, compressed_ranges, padding_amt, output_ref);
		std::cerr << "expressing variants with " << compressed_ranges.size() << " sequences." << std::endl;
	}
}
