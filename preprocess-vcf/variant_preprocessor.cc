/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/transform.hpp>
#include <vcf2multialign/check_overlapping_non_nested_variants.hh>
#include <vcf2multialign/preprocess/variant_preprocessor.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void variant_preprocessor::output_sample_names()
	{
		// XXX I don’t remember why vcf_reader outputs the sample names as a map with names as keys and positions as indices but it should be fixed unless a good reason is found not to.
		std::vector <std::string> sample_names;
		for (auto const & [sample_name, idx1] : m_reader->sample_names())
			sample_names[idx1 - 1] = sample_name; // Copy.
	
		for (auto const &sample_name : sample_names)
			*m_ostream << '\t' << sample_name;
	}
	
	
	void variant_preprocessor::output_headers()
	{
		*m_ostream << "##fileformat=VCFv4.3\n";
		*m_ostream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		for (auto const &[key, contig] : m_reader->metadata().contig())
		{
			contig.output_vcf(*m_ostream);
			*m_ostream << '\n';
		}
		
		*m_ostream << "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		output_sample_names();
		*m_ostream << '\n';
	}
	
	
	void variant_preprocessor::process_and_output(std::vector <std::string> const &field_names_for_filter_by_assigned)
	{
		output_headers();
		
		m_reader->reset();
		m_reader->set_parsed_fields(lb::vcf_field::ALL);
		m_overlap_start = 0;
		m_overlap_end = 0;
		
		// Get the field descriptors needed for accessing the values.
		libbio::vcf_info_field_end const	*end_field{};
		m_reader->get_info_field_ptr("END",	end_field);
		
		// Determine the fields used for filtering.
		std::vector <lb::vcf_info_field_base *> filter_by_assigned;
		{
			auto const &fields(m_reader->info_fields());
			for (auto const &name : field_names_for_filter_by_assigned)
			{
				auto const it(fields.find(name));
				if (fields.end() == it)
				{
					std::cerr << "WARNING: Did not find a field for identifier “" << name << "”.\n";
					continue;
				}
				
				filter_by_assigned.emplace_back(it->second.get());
			}
		}
		
		// Process the variants.
		bool should_continue(false);
		do {
			m_reader->fill_buffer();
			should_continue = m_reader->parse(
				[
					this,
					end_field,
					&filter_by_assigned
				](lb::transient_variant const &var) -> bool
				{
					// Check the chromosome name.
					if (var.chrom_id() != m_chromosome_name)
						return true;
					
					auto const lineno(var.lineno());
					auto const var_pos(var.zero_based_pos());
					libbio_always_assert_msg(var_pos < m_reference->size(), "Variant position on line ", lineno, " greater than reference length (", m_reference->size(), ").");
					
					for (auto const &alt : var.alts())
					{
						// FIXME: log if skipped.
						if (!can_handle_variant_alt(alt))
							return true;
					}
					
					// Filter.
					for (auto const *field_ptr : filter_by_assigned)
					{
						// FIXME: log if skipped.
						if (field_ptr->has_value(var))
							return true;
					}
					
					// Variant passes the checks, handle it.
					{
						auto const var_pos(var.zero_based_pos());
						
						// As long as there are parallel paths, add to the overlap stack.
						if (var_pos < m_overlap_end + m_minimum_subgraph_distance)
							m_seen_parallel_paths = true;
						else if (m_seen_parallel_paths) // true == m_overlap_end + m_minimum_subgraph_distance <= var_pos
						{
							// Reset.
							m_seen_parallel_paths = false;
							process_overlap_stack();
							m_overlap_start = var_pos;
						}
						
						m_overlapping_variants.emplace_back(var);
						auto const var_end(lb::variant_end_pos(var, *end_field));
						m_overlap_end = std::max(m_overlap_end, var_end);
					}
					return true;
				}
			);
		} while (should_continue);
		
		// Process the remaining variants.
		process_overlap_stack();
	}
	
	
	// Retrieve the accumulated group of variants and pass them to the worker thread for processing.
	void variant_preprocessor::process_overlap_stack()
	{
		// Fast path: empty stack.
		if (m_overlapping_variants.empty())
			return;
		
		// Fast path: single variant.
		if (1 == m_overlapping_variants.size())
		{
			output_single_variant(m_overlapping_variants.front());
			m_overlapping_variants.clear();
			return;
		}
		
		// Slow path: overlapping variants.
		m_sample_sorter.prepare_for_next_subgraph();
		for (auto const &var : m_overlapping_variants)
		{
			for (std::size_t i(0), count(var.alts().size()); i < count; ++i)
				m_sample_sorter.sort_by_variant_and_alt(var, i);
		}
	
		path_sorted_variant psv;
	
		auto const path_count(m_sample_sorter.path_count());
		psv.set_paths_by_sample(m_sample_sorter.paths_by_sample());
		psv.reserve_memory_for_paths(path_count);
		psv.invert_paths_by_sample();
	
		// Determine the ALT sequence for each path by taking a representative for each path.
		lb::packed_matrix <1> alt_sequences_by_path(m_overlapping_variants.size(), path_count);
		{
			// Iterate over the variants and ALTs.
			std::size_t row_idx(0);
			for (auto const &var : m_overlapping_variants)
			{
				// Iterate over the sample representatives.
				auto const *gt_field(get_variant_format(var).gt);
				std::size_t path_idx(0);
				for (auto const &sample_numbers : psv.samples_by_path())
				{
					libbio_always_assert(!sample_numbers.empty());
					auto const sample_idx(sample_numbers.front());
					auto const [donor_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(sample_idx);
				
					// Store the GT value.
					auto const alt((*gt_field)(var.samples()[donor_idx])[chr_idx].alt);
					alt_sequences_by_path(row_idx, path_idx).fetch_or(alt);
				
					++path_idx;
				}
			
				++row_idx;
			}
		}
	
		output_sorted(psv, alt_sequences_by_path);
		m_overlapping_variants.clear();
	}
	
	
	void variant_preprocessor::output_sorted(path_sorted_variant const &psv, lb::packed_matrix <1> const &alt_sequences_by_path)
	{
		// Generate the REF and ALT values for outputting one variant.
		// The path numbers may be used for the GT values.
		std::string_view reference_sv(m_reference->data(), m_reference->size());
		std::string ref_string;
		std::vector <std::string> alt_strings_by_path(alt_sequences_by_path.number_of_columns());
		auto const path_count(m_sample_sorter.path_count());
	
		{
			// Iterate over the variants and ALTs.
			std::size_t ref_start(m_overlap_start);
			std::size_t row_idx(0);
			for (auto const &var : m_overlapping_variants)
			{
				for (std::size_t path_idx(0); path_idx < path_count; ++path_idx)
				{
					// Output the part of the reference that is between the variants.
					alt_strings_by_path[path_idx] += reference_sv.substr(ref_start, var.zero_based_pos() - ref_start);
				
					auto const alt_idx(alt_sequences_by_path(row_idx, path_idx));
					if (0 == alt_idx)
					{
						// Output REF.
						alt_strings_by_path[path_idx] += reference_sv.substr(var.zero_based_pos(), lb::variant_end_pos(var, *m_end_field));
					}
					else
					{
						// Output based on the ALT type.
						auto const &alt(var.alts()[alt_idx]);
						auto const svt(alt.alt_sv_type);
						switch (svt)
						{
							case lb::sv_type::NONE:
								// Output the ALT sequence.
								alt_strings_by_path[path_idx] += alt.alt;
								break;
								
							case lb::sv_type::DEL:
							case lb::sv_type::DEL_ME:
								// Don’t output anything.
								break;
								
							case lb::sv_type::INS:
							case lb::sv_type::DUP:
							case lb::sv_type::INV:
							case lb::sv_type::CNV:
							case lb::sv_type::DUP_TANDEM:
							case lb::sv_type::INS_ME:
							case lb::sv_type::UNKNOWN:
							default:
								throw std::runtime_error("Unexpected ALT type");
								break;
						}
					}
				}
			
				ref_string += var.ref();
				ref_start += var.ref().size();
				++row_idx;
			}
		}
		
		// Check for empty ALT strings, substitute with <DEL>.
		for (auto &alt : alt_strings_by_path)
		{
			if (alt.empty())
				alt = "<DEL>";
		}
		
		output_combined_variants(psv, ref_string, alt_strings_by_path);
	}
	
	
	void variant_preprocessor::output_single_variant(lb::variant const &var) const
	{
		// #CHROM, POS, ID, REF
		*m_ostream
			<< m_chromosome_name << '\t'
			<< var.pos() << '\t'
			<< ".\t"
			<< var.ref() << '\t';
		
		// ALT
		ranges::copy(
			var.alts() | ranges::view::transform([](auto const &alt) -> std::string const & { return alt.alt; }),
			ranges::make_ostream_joiner(*m_ostream, "\t")
		);
		*m_ostream << '\t';
		
		// QUAL, FILTER, INFO, FORMAT
		*m_ostream
			<< ".\t"
			<< "PASS\t"
			<< ".\t"
			<< "GT";
		
		// Samples.
		auto const *gt_field(get_variant_format(var).gt);
		for (auto const &sample : var.samples())
		{
			*m_ostream << '\t';
			gt_field->output_vcf_value(*m_ostream, sample);
		}
		*m_ostream << '\n';
	}
	
	
	void variant_preprocessor::output_combined_variants(
		path_sorted_variant const &psv,
		std::string const &ref_string,
		std::vector <std::string> const &alt_strings_by_path
	)
	{
		// Output combined variants as a VCF record.
		// Our sample order should match that of the VCF file.
		// Copy the values to make sure, though.
		auto const &paths_by_sample(psv.paths_by_sample());
		std::fill(m_output_gt.begin(), m_output_gt.end(), 0);
		for (std::size_t i(0), count(paths_by_sample.size()); i < count; ++i)
		{
			auto const [chr_idx, sample_idx] = m_sample_indexer.donor_and_chr_idx(i);
			m_output_gt(chr_idx, sample_idx) = paths_by_sample[i];
		}
		
		// #CHROM, POS, ID, REF
		*m_ostream
			<< m_chromosome_name << '\t'
			<< m_overlap_start << '\t'
			<< ".\t"
			<< ref_string << '\t';
		
		// ALT
		ranges::copy(alt_strings_by_path, ranges::make_ostream_joiner(*m_ostream, "\t"));
		*m_ostream << '\t';
		
		// QUAL, FILTER, INFO, FORMAT
		*m_ostream
			<< ".\t"
			<< "PASS\t"
			<< ".\t"
			<< "GT";
		
		// FIXME: add row and column ranges to libbio::matrix and libbio::packed_matrix.
		for (std::size_t i(0), count(m_output_gt.number_of_columns()); i < count; ++i)
		{
			*m_ostream << '\t';
			
			auto const &column(m_output_gt.const_column(i));
			ranges::copy(column, ranges::make_ostream_joiner(*m_ostream, "|"));
		}
		*m_ostream << '\n';
	}
}
