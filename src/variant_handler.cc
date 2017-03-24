/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/util.hh>
#include <vcf2multialign/variant_handler.hh>


namespace vcf2multialign {
	
	// Fill the streams with '-'.
	void variant_handler::fill_streams(haplotype_ptr_map &haplotypes, size_t const fill_amt) const
	{
		for (auto &kv : haplotypes)
		{
			for (auto h_ptr : kv.second)
			{
				if (h_ptr)
				{
					std::ostream_iterator <char> it(h_ptr->output_stream);
					std::fill_n(it, fill_amt, '-');
				}
			}
		}
	}

	
	// Fill the streams with reference.
	void variant_handler::output_reference(std::size_t const output_start_pos, std::size_t const output_end_pos)
	{
		if (output_start_pos == output_end_pos)
			return;
		
		always_assert(output_start_pos < output_end_pos, "Bad offset order");

		char const *ref_begin(m_reference->data());
		auto const output_start(ref_begin + output_start_pos);
		for (auto &kv : m_ref_haplotype_ptrs)
		{
			for (auto h_ptr : kv.second)
			{
				if (h_ptr)
				{
					auto &h(*h_ptr);
					always_assert(h.current_pos == output_start_pos, "Unexpected position");
					
					h.output_stream.write(output_start, output_end_pos - output_start_pos);
					h.current_pos = output_end_pos;
				}
			}
		}
	}

	
	std::size_t variant_handler::process_overlap_stack(size_t const var_pos)
	{
		std::size_t retval(0);
		while (true)
		{
			// NOTE: for libc++, use p overlap_stack.c in the debugger (not p overlap_stack).
			variant_overlap &vo(m_overlap_stack.top());
			if (var_pos < vo.end_pos)
				break;

			// The sequence was output up to vo's start_pos when it was added to the stack.
			// Check the current position and output from there.
			auto const output_start_pos(vo.current_pos);
			auto const output_end_pos(vo.end_pos);

			// Return the end position of the last variant overlap.
			retval = output_end_pos;
			
			// Output reference from 5' direction up to vo.end_pos.
			output_reference(vo.current_pos, vo.end_pos);
			
			// Add the amount output to the heaviest path length.
			vo.heaviest_path_length += vo.end_pos - vo.current_pos;

			// Also output ALTs. At the same time,
			// compare the heaviest path length to the lengths of every alternative.
			auto heaviest_path_length(vo.heaviest_path_length);
			for (auto &kv : vo.alt_haplotypes)
			{
				auto const &alt(kv.first);
				auto &haplotypes(kv.second);
				for (auto &kv : haplotypes)
				{
					for (auto h_ptr : kv.second)
					{
						if (h_ptr)
						{
							auto &h(*h_ptr);
							always_assert(h.current_pos <= output_start_pos, "Unexpected position");
							
							h.output_stream << alt;
							h.current_pos = output_end_pos;
						}
					}
				}

				heaviest_path_length = std::max(heaviest_path_length, alt.size());
			}

			// Update current_pos to match the position up to which the sequence was output.
			vo.current_pos = output_end_pos;

			// Fill the shorter sequences. vo.heaviest_path_length considers REF only.
			if (vo.heaviest_path_length < heaviest_path_length)
			{
				auto const fill_amt(heaviest_path_length - vo.heaviest_path_length);
				fill_streams(m_ref_haplotype_ptrs, fill_amt);
			}

			// Also fill the shorter alternatives.
			for (auto &kv : vo.alt_haplotypes)
			{
				auto const &alt(kv.first);
				auto const alt_length(alt.size());
				if (alt_length < heaviest_path_length)
				{
					auto &haplotypes(kv.second);
					auto const fill_amt(heaviest_path_length - alt_length);
					fill_streams(haplotypes, fill_amt);
				}
			}

			// Since the variants were written to the haplotypes, they can now be updated with the reference again.
			for (auto &kv : vo.alt_haplotypes)
			{
				auto &haplotypes(kv.second);
				for (auto &kv : haplotypes)
				{
					auto const &sample_name(kv.first);
					auto &alt_ptrs(kv.second);
					auto &ref_ptrs(m_ref_haplotype_ptrs[sample_name]);
					
					for (size_t i(0), count(alt_ptrs.size()); i < count; ++i)
					{
						always_assert(! (alt_ptrs[i] && ref_ptrs[i]), "Inconsistent haplotype pointers");
						
						// Use ADL.
						using std::swap;
						if (alt_ptrs[i])
							swap(alt_ptrs[i], ref_ptrs[i]);
					}
				}
			}

			// Update current_pos and heaviest path length with the result.
			if (1 < m_overlap_stack.size())
			{
				m_overlap_stack.pop();
				auto &previous_overlap(m_overlap_stack.top());
				previous_overlap.current_pos = output_end_pos;
				previous_overlap.heaviest_path_length += heaviest_path_length;
			}
			else
			{
				// If the handled variant_overlap was the last one, exit the loop.
				break;
			}
		}

		return retval;
	}

		
	void variant_handler::handle_variant(variant &var)
	{
		m_skipped_samples.clear();
		m_counts_by_alt.clear();
		m_non_ref_totals.reset();
		
		auto const lineno(var.lineno());
		if (0 != m_skipped_variants->count(lineno))
			return;
		
		// Preprocess the ALT field to check that it can be handled.
		{
			// Check that the alt sequence is something that can be handled.
			m_valid_alts.clear();
			std::size_t i(0);
			for (auto const &alt : var.alts())
			{
				++i;
				for (auto const c : alt)
				{
					if (! ('A' == c || 'C' == c || 'G' == c || 'T' == c || 'N' == c))
					{
						std::cerr << "Unable to handle ALT field " << alt << " on line " << lineno << '.' << std::endl;
						goto loop_end;
					}
				}
				
				m_valid_alts.emplace(i);
				
			loop_end:
				;
			}
			
			if (m_valid_alts.empty())
				return;
		}
		
		size_t const var_pos(var.zero_based_pos());
		
		auto const var_ref(var.ref());
		auto const var_ref_size(var_ref.size());
		auto const var_alts(var.alts());
		
		// If var is beyond previous_variant.end_pos, handle the variants on the stack
		// until a containing variant is found or the bottom of the stack is reached.
		process_overlap_stack(var_pos);
		
		// Use the previous variant's range to determine the output sequence.
		auto &previous_variant(m_overlap_stack.top());
		
		// Check that if var is before previous_variant.end_pos, it is also completely inside it.
		always_assert(
			(previous_variant.end_pos <= var_pos) || (var_pos + var_ref_size <= previous_variant.end_pos),
			[lineno, &previous_variant, var_pos, var_ref_size](){
				std::cerr
				<< "Invalid variant inclusion."
				<< "\n\tlineno:\t" << lineno
				<< "\n\tprevious_variant.lineno:\t" << previous_variant.lineno
				<< "\n\tvar_pos:\t" << var_pos
				<< "\n\tvar_ref_size:\t" << var_ref_size
				<< "\n\tprevious_variant.start_pos:\t" << previous_variant.start_pos
				<< "\n\tprevious_variant.end_pos:\t" << previous_variant.end_pos
				<< std::endl;
			}
		);
		
		// Output reference from 5' direction up to var_pos.
		output_reference(previous_variant.current_pos, var_pos);
		
		// Add the amount output to the heaviest path length.
		previous_variant.heaviest_path_length += var_pos - previous_variant.current_pos;
		
		// Update current_pos to match the position up to which the sequence was output.
		// Also add the length to the heaviest path length.
		previous_variant.current_pos = var_pos;
		
		// Find haplotypes that have the variant.
		// First make sure that all valid alts are listed in m_alt_haplotypes.
		for (auto const alt_idx : m_valid_alts)
		{
			auto const &alt_str(var_alts[alt_idx - 1]);
			m_alt_haplotypes[alt_str];
		}
		m_alt_haplotypes[*m_null_allele_seq];
				
		for (auto const &kv : *m_all_haplotypes)
		{
			auto const sample_no(kv.first);

			// Get the sample.
			auto const sample(var.sample(sample_no));
			
			// Handle the genotype.
			uint8_t chr_idx(0);
			for (auto const gt : sample.get_genotype())
			{
				auto const alt_idx(gt.alt);
				auto const is_phased(gt.is_phased);
				always_assert(0 == chr_idx || is_phased, "Variant file not phased");

				if (0 != alt_idx && 0 != m_valid_alts.count(alt_idx))
				{
					auto &ref_ptrs(m_ref_haplotype_ptrs.find(sample_no)->second); // Has nodes for every sample_no.

					std::string const *alt_ptr{m_null_allele_seq};
					if (NULL_ALLELE != alt_idx)
						alt_ptr = &var_alts[alt_idx - 1];
					
					haplotype_ptr_map &alt_ptrs_by_sample(m_alt_haplotypes[*alt_ptr]);
					auto it(alt_ptrs_by_sample.find(sample_no));
					if (alt_ptrs_by_sample.end() == it)
					{
						it = alt_ptrs_by_sample.emplace(
							std::piecewise_construct,
							std::forward_as_tuple(sample_no),
							std::forward_as_tuple(ref_ptrs.size(), nullptr)
						).first;
					}
					auto &alt_ptrs(it->second);
					
					if (ref_ptrs[chr_idx])
					{
						// Use ADL.
						using std::swap;
						swap(alt_ptrs[chr_idx], ref_ptrs[chr_idx]);
						++m_non_ref_totals.handled_count;
						++m_counts_by_alt[alt_idx].handled_count;
					}
					else
					{
						if (m_overlapping_alts.insert(lineno).second)
						{
							std::cerr << "Overlapping alternatives on line " << lineno
							<< " for sample " << sample_no << ':' << (int) chr_idx
							<< " (and possibly others); skipping when needed." << std::endl;
						}
						
						if (m_error_logger->is_logging_errors())
							m_skipped_samples.emplace_back(sample_no, alt_idx, chr_idx);
					}
					
					++m_non_ref_totals.total_count;
					++m_counts_by_alt[alt_idx].total_count;
				}
				++chr_idx;
			}
		}
		
		// Report errors if needed.
		if (m_error_logger->is_logging_errors())
		{
			for (auto const &s : m_skipped_samples)
				m_error_logger->log_overlapping_alternative(lineno, s.sample_no, s.chr_idx, m_counts_by_alt[s.alt_idx], m_non_ref_totals);
		}
		
		// Create a new variant_overlap.
		auto const var_end(var_pos + var_ref_size);
		auto const previous_end_pos(previous_variant.end_pos);
		variant_overlap overlap(var_pos, var_pos, var_end, 0, lineno, m_alt_haplotypes);
		if (var_pos < previous_end_pos)
		{
			// Add the current variant to the stack.
			m_overlap_stack.emplace(std::move(overlap));
		}
		else
		{
			// Replace the top of the stack with the current variant.
			// Use ADL.
			using std::swap;
			swap(m_overlap_stack.top(), overlap);
		}
		
		++m_i;
		if (0 == m_i % 50000)
			std::cerr << "Handled " << m_i << " variants…" << std::endl;
	}
	
	
	void variant_handler::finish()
	{
		// Fill the remaining part with reference.
		m_error_logger->flush();
		std::cerr << "Filling with the reference…" << std::endl;
		auto const ref_size(m_reference->size());
		auto const output_end_pos(process_overlap_stack(ref_size));
		
		char const *ref_begin(m_reference->data());
		for (auto &kv : *m_all_haplotypes)
		{
			for (auto &h : kv.second)
			{
				auto const output_len(ref_size - h.current_pos);
				h.output_stream.write(ref_begin + h.current_pos, output_len);
			}
		}
		
		m_finish_callback();
	}
	
	
	void variant_handler::process_variants(haplotype_map &all_haplotypes)
	{
		while (!m_overlap_stack.empty())
			m_overlap_stack.pop();
		
		m_ref_haplotype_ptrs.clear();
		m_overlap_stack.emplace(0, 0, 0, 0, 0);
		m_all_haplotypes = &all_haplotypes;
		m_i = 0;
		
		// All haplotypes initially have the reference sequence.
		for (auto &kv : *m_all_haplotypes)
		{
			auto const sample_no(kv.first);
			auto &haplotype_vector(kv.second);
			auto &haplotype_ptr_vector(m_ref_haplotype_ptrs[sample_no]);
			
			auto const count(haplotype_vector.size());
			haplotype_ptr_vector.resize(count);
			for (size_t i(0); i < count; ++i)
				haplotype_ptr_vector[i] = &haplotype_vector[i];
		}
		
		auto &reader(m_variant_buffer.reader());
		reader.reset();
		reader.set_parsed_fields(vcf_field::ALL);
		
		dispatch_async_f <decltype(m_variant_buffer), &variant_buffer::read_input>(*m_parsing_queue, &m_variant_buffer);

	}
}
