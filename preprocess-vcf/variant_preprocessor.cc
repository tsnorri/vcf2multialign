/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/transform.hpp>
#include <vcf2multialign/can_handle_variant_alts.hh>
#include <vcf2multialign/preprocess/variant_preprocessor.hh>
#include <vcf2multialign/variant_format.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	void variant_preprocessor::update_sample_names()
	{
		// XXX I don’t remember why vcf_reader outputs the sample names as a map with names as keys and positions as indices but it should be fixed unless a good reason is found not to.
		auto const &sample_name_map(m_reader->sample_names());
		m_sample_names.resize(sample_name_map.size());
		for (auto const & [sample_name, idx1] : sample_name_map)
			m_sample_names[idx1 - 1] = sample_name; // Copy.
	}
	
	
	void variant_preprocessor::output_headers()
	{
		*m_ostream << "##fileformat=VCFv4.3\n";
		*m_ostream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		*m_ostream << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
		m_output_lineno = 3;
		for (auto const &[key, contig] : m_reader->metadata().contig())
		{
			contig.output_vcf(*m_ostream);
			++m_output_lineno;
		}
		
		*m_ostream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		for (auto const &sample_name : m_sample_names)
			*m_ostream << '\t' << sample_name;
		*m_ostream << '\n';
		++m_output_lineno;
	}
	
	
	void variant_preprocessor::process_and_output(std::vector <std::string> const &field_names_for_filter_by_assigned)
	{
		update_sample_names();
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
		std::size_t processed_count(0);
		do {
			m_reader->fill_buffer();
			should_continue = m_reader->parse(
				[
					this,
					end_field,
					&filter_by_assigned,
				 	&processed_count
				](lb::transient_variant const &var) -> bool
				{
					auto const lineno(var.lineno());
					auto const var_pos(var.zero_based_pos());

					// Check the chromosome name.
					if (var.chrom_id() != m_chromosome_name)
						goto end;
					
					libbio_always_assert_msg(var_pos < m_reference->size(), "Variant position on line ", lineno, " greater than reference length (", m_reference->size(), ").");

					if (!can_handle_variant_alts(var))
					{
						std::cerr << "Line " << lineno << ": Variant has no ALTs that could be handled.\n";
						goto end;
					}
					
					// Filter.
					for (auto const *field_ptr : filter_by_assigned)
					{
						if (field_ptr->has_value(var))
						{
							std::cerr << "Line " << lineno << ": Variant has the field '" << field_ptr->get_metadata()->get_id() << "' set; skipping.\n";
							goto end;
						}
					}
					
					// Check the reference seuqence.
					{
						auto const &ref_col(var.ref());
						std::string_view const ref_sub(m_reference->data() + var_pos, ref_col.size());
						if (ref_col != ref_sub)
						{
							std::cerr << "WARNING: reference column mismatch on line " << lineno << ": expected '" << ref_sub << "', got '" << ref_col << "'\n";
							goto end;
						}
					}
					
					// Variant passes the checks, handle it.
					{
						if (m_overlap_end + m_minimum_subgraph_distance <= var_pos)
						{
							process_overlap_stack();
							m_overlap_start = var_pos;
						}
						
						m_overlapping_variants.emplace_back(var);
						auto const var_end(lb::variant_end_pos(var, *end_field));
						m_overlap_end = std::max(m_overlap_end, var_end);
					}
					
				end:
					++processed_count;
					if (0 == processed_count % 100000)
						std::cerr << "Handled " << processed_count << " variants…\n";
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
				m_sample_sorter.sort_by_variant_and_alt(var, 1 + i);
		}
	
		path_sorted_variant psv;
	
		auto const path_count(m_sample_sorter.path_count());
		std::cerr << "Line " << m_overlapping_variants.front().lineno() << ": Overlapping variants: " << m_overlapping_variants.size() << " path count: " << path_count << '\n';
		psv.set_paths_by_sample(m_sample_sorter.paths_by_sample());
		psv.reserve_memory_for_paths(path_count);
		psv.invert_paths_by_sample();
	
		// Determine the ALT sequence for each path by taking a representative for each path.
		// FIXME: add a check (other than the assertion in the matrix) for the ALT value’s limit.
		lb::packed_matrix <8> alt_sequences_by_path(m_overlapping_variants.size(), path_count);
		{
			// Iterate over the variants and ALTs.
			std::size_t row_idx(0);
			for (auto const &var : m_overlapping_variants)
			{
				try
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
						auto const &sample(var.samples()[donor_idx]);
						auto const &genotypes((*gt_field)(sample));
						auto const &genotype(genotypes[chr_idx]);
						auto const alt(genotype.alt);
						if (lb::sample_genotype::NULL_ALLELE == alt)
							std::cerr << "Line " << var.lineno() << ": Ignoring null allele in sample " << donor_idx << ':' << chr_idx << '\n';
						else
							alt_sequences_by_path(row_idx, path_idx).fetch_or(alt);
						//std::cerr << "4 row_idx: " << row_idx << " path_idx: " << path_idx << " alt: " << alt << '\n';
					
						++path_idx;
					}
				
					++row_idx;
				}
				catch (std::exception const &exc)
				{
					std::cerr << "FATAL: Caught an exception while handling the following variant:" << std::endl;
					lb::output_vcf(std::cerr, var);
					std::cerr << std::endl;
					throw exc;
				}
			}
		}
	
		output_sorted(psv, alt_sequences_by_path);
		m_overlapping_variants.clear();
	}
	
	
	void variant_preprocessor::output_sorted(path_sorted_variant const &psv, lb::packed_matrix <8> const &alt_sequences_by_path)
	{
		// Generate the REF and ALT values for outputting one variant.
		// The path numbers may be used for the GT values.
		std::string_view reference_sv(m_reference->data(), m_reference->size());
		std::vector <std::string> alt_strings_by_path(alt_sequences_by_path.number_of_columns());
		std::vector <std::size_t> alt_string_output_positions(alt_sequences_by_path.number_of_columns(), m_overlap_start);
		
		// Iterate over the variants and ALTs.
		std::size_t row_idx(0);
		for (auto const &pair : m_overlapping_variants | ranges::view::sliding(2))
		{
			auto const &var(pair[0]);
			auto const &next_var(pair[1]);
			auto const next_var_pos(next_var.zero_based_pos());
			
			output_sorted_handle_paths(
				row_idx,
				var,
				next_var_pos,
				reference_sv,
				alt_sequences_by_path,
				psv.samples_by_path(),
				alt_strings_by_path,
				alt_string_output_positions
			);
			
			++row_idx;
		}
		
		// Handle the last variant.
		{
			auto const &var(m_overlapping_variants.back());
			auto const end_pos(lb::variant_end_pos(var, *m_end_field));
			output_sorted_handle_paths(
				row_idx,
				var,
				end_pos,
				reference_sv,
				alt_sequences_by_path,
				psv.samples_by_path(),
				alt_strings_by_path,
				alt_string_output_positions
			);
			
			// Fill up to the overlap end.
			for (auto const pair : ranges::view::zip(alt_string_output_positions, alt_strings_by_path))
			{
				auto const &output_pos(std::get <0>(pair));
				auto &alt_str(std::get <1>(pair));
				auto const ref_part(reference_sv.substr(output_pos, end_pos - output_pos));
				alt_str += ref_part;
			}
		}
		
		// Check for empty ALT strings, substitute with <DEL>.
		for (auto &alt : alt_strings_by_path)
		{
			if (alt.empty())
				alt = "<DEL>";
		}
		
		std::cerr << "Writing to line " << m_output_lineno << '\n';
		output_combined_variants(psv, alt_strings_by_path);
	}
	
	
	void variant_preprocessor::output_sorted_handle_paths(
		std::size_t const row_idx,
		libbio::variant const &var,
		std::size_t const next_var_pos,
		std::string_view const &reference_sv,
		lb::packed_matrix <8> const &alt_sequences_by_path,
		sample_map const &samples_by_path,
		std::vector <std::string> &alt_strings_by_path,
		std::vector <std::size_t> &alt_string_output_positions
	) const
	{
		auto const var_pos(var.zero_based_pos());
		auto const end_pos(lb::variant_end_pos(var, *m_end_field));
		
		auto const path_count(m_sample_sorter.path_count());
		for (std::size_t path_idx(0); path_idx < path_count; ++path_idx)
		{
			auto const alt_end(alt_string_output_positions[path_idx]);
			auto const alt_idx(alt_sequences_by_path(row_idx, path_idx));
			// Skip the variants that would overlap with the previous variant on the path in question.
			if (var_pos < alt_end)
			{
				if (0 != alt_idx)
				{
					std::cerr << "Line " << var.lineno() << ": Overlap on path " << path_idx << " in sample(s) ";
					bool is_first(true);
					for (auto const sample_no : samples_by_path[path_idx])
					{
						if (!is_first)
							std::cerr << ", ";
						auto const [donor_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(sample_no);
						std::cerr << m_sample_names[donor_idx] << ':' << +chr_idx << " (" << sample_no << ')';
						is_first = false;
					}
					std::cerr << "; skipping.\n";
				}
				continue;
			}
			
			// Output the part of the reference that is between the variants.
			// This corresponds to joining paths in the DAG.
			alt_strings_by_path[path_idx] += reference_sv.substr(alt_end, var_pos - alt_end);
			
			libbio_assert(0 != path_idx || 0 == alt_idx); // path_idx 0 is supposed to have the reference sequence.
			if (0 == alt_idx)
			{
				// Output REF.
				alt_strings_by_path[path_idx] += reference_sv.substr(var_pos, next_var_pos - var_pos);
				alt_string_output_positions[path_idx] = next_var_pos;
			}
			else
			{
				// Output based on the ALT type.
				auto const &alt(var.alts()[alt_idx - 1]);
				auto const svt(alt.alt_sv_type);
				switch (svt)
				{
					case lb::sv_type::NONE:
						// Output the ALT sequence.
						alt_strings_by_path[path_idx] += alt.alt;
						alt_string_output_positions[path_idx] = end_pos;
						break;
						
					case lb::sv_type::DEL:
					case lb::sv_type::DEL_ME:
						// Don’t output anything.
						alt_string_output_positions[path_idx] = end_pos;
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
	}
	
	
	void variant_preprocessor::output_single_variant(lb::variant const &var)
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
		++m_output_lineno;
	}
	
	
	void variant_preprocessor::output_combined_variants(
		path_sorted_variant const &psv,
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
			auto const [sample_idx, chr_idx] = m_sample_indexer.donor_and_chr_idx(i);
			m_output_gt(chr_idx, sample_idx) = paths_by_sample[i];
		}
		
		// #CHROM, POS, ID, REF
		*m_ostream
			<< m_chromosome_name << '\t'
			<< 1 + m_overlap_start << '\t'	// m_overlap_start is zero-based.
			<< ".\t"
			<< alt_strings_by_path.front() << '\t';
		
		// ALT
		ranges::copy(alt_strings_by_path | ranges::view::drop(1), ranges::make_ostream_joiner(*m_ostream, ","));
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
			ranges::copy(
				column | ranges::view::transform([](auto const &gt_idx){ return +gt_idx; }),
				ranges::make_ostream_joiner(*m_ostream, "|")
			);
		}
		*m_ostream << '\n';
		++m_output_lineno;
	}
}
