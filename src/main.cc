/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
*/

#include <boost/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <cmath>
#include <fcntl.h>
#include <iostream>
#include <map>
#include <Variant.h>
#include <vcf2multialign/cxx14compat.hh>
#include <vcf2multialign/fasta_reader.hh>
#include <set>
#include <stack>
#include "cmdline.h"

namespace ios	= boost::iostreams;
namespace v2m	= vcf2multialign;
namespace vl	= vcflib;


namespace {
	typedef std::vector <char> vector_type;

	typedef ios::stream <ios::file_descriptor_source>	file_istream;
	typedef ios::stream <ios::file_descriptor_sink>		file_ostream;

	typedef std::set <std::string> variant_set;
	typedef std::set <std::string> sample_set;
	
	typedef std::map <std::string, size_t> ploidy_map;
	

	struct haplotype
	{
		size_t current_pos{0};
		file_ostream output_stream;
	};
	
	
	typedef std::map <
		std::string,				// Sample name
		std::vector <haplotype>		// All haplotype sequences
	> haplotype_map;
	
	
	typedef std::map <
		std::string,				// Sample name
		std::vector <haplotype *>	// Haplotype sequences by chromosome index
	> haplotype_ptr_map;
	
	
	typedef std::map <
		std::string,				// ALT
		haplotype_ptr_map
	> alt_map;
	
	
	struct variant_overlap
	{
		size_t					start_pos{0};
		size_t					current_pos{0};
		size_t					end_pos{0};
		size_t					heaviest_path_length{0};
		alt_map					alt_haplotypes;
		
		variant_overlap(size_t const s, size_t const c, size_t const e, size_t const h):
			start_pos(s),
			current_pos(c),
			end_pos(e),
			heaviest_path_length(h)
		{
			if (! (start_pos <= end_pos))
				throw std::runtime_error("Bad offset order");
		}
		
		variant_overlap(size_t const s, size_t const c, size_t const e, size_t const h, alt_map &alts):
			start_pos(s),
			current_pos(c),
			end_pos(e),
			heaviest_path_length(h),
			alt_haplotypes(std::move(alts))
		{
			if (! (start_pos <= end_pos))
				throw std::runtime_error("Bad offset order");
		}
	};
	
	
	typedef std::stack <variant_overlap> overlap_stack_type;
	
	
	void handle_file_error(char const *fname)
	{
		char const *errmsg(strerror(errno));
		std::cerr << "Got an error while trying to open '" << fname << "': " << errmsg << std::endl;
		exit(EXIT_FAILURE);
	}


	void open_file_for_reading(char const *fname, file_istream &stream)
	{
		int fd(open(fname, O_RDONLY));
		if (-1 == fd)
			handle_file_error(fname);

		ios::file_descriptor_source source(fd, ios::close_handle);
		stream.open(source);
	}


	void open_file_for_writing(char const *fname, file_ostream &stream, bool const should_overwrite)
	{
		int fd(0);
		if (should_overwrite)
			fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		else
			fd = open(fname, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);

		if (-1 == fd)
			handle_file_error(fname);

		ios::file_descriptor_sink sink(fd, ios::close_handle);
		stream.open(sink);
	}


	void open_files_for_writing(
		sample_set const &sample_names,
		haplotype_map &haplotypes,
		ploidy_map const &ploidy,
		bool const should_overwrite
	)
	{
		haplotypes.clear(); // Closes files since ios::close_handle was given.
		for (auto const &sample_name : sample_names)
		{
			auto const current_ploidy(ploidy.find(sample_name)->second);
			
			// Since file_ostream is not movable, check first if the vector has already been created.
			// If not, create it with the exact size to avoid resizing and thus moving later.
			auto it(haplotypes.find(sample_name));
			if (haplotypes.end() == it)
			{
				it = haplotypes.emplace(
					std::piecewise_construct,
					std::forward_as_tuple(sample_name),
					std::forward_as_tuple(current_ploidy)
				).first;
			}
			
			auto &haplotype_vec(it->second);
			
			for (size_t i(1); i <= current_ploidy; ++i)
			{
				auto const fname(boost::str(boost::format("%s-%u") % sample_name % i));
				open_file_for_writing(fname.c_str(), haplotype_vec[i - 1].output_stream, should_overwrite);
			}
		}
	}

	
	size_t zero_based_position_safe(vl::Variant const &var)
	{
		auto const temp_pos(var.zeroBasedPosition());
		if (temp_pos < 0)
			throw std::runtime_error("Unexpected position");
		
		return static_cast <size_t>(temp_pos);
	}
	
	
	void fill_streams(haplotype_ptr_map &haplotypes, size_t const fill_amt)
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
	
	
	// For using std::stringstream as output streams.
#if 0
	extern void print_seq(haplotype_map const &all_haplotypes)
	{
		for (auto &kv : all_haplotypes)
		{
			auto const &sample(kv.first);
			std::size_t i(1);
			for (auto &h : kv.second)
			{
				std::cerr << sample << '-' << i << ": " << h.output_stream.str() << std::endl;
				++i;
			}
		}
		std::cerr << std::endl;
	}
#endif


	void output_reference(
		haplotype_ptr_map &ref_haplotype_ptrs,
		vector_type const &reference,
		std::size_t const output_start_pos,
		std::size_t const output_end_pos
	)
	{
		if (output_start_pos == output_end_pos)
			return;
		
		if (! (output_start_pos < output_end_pos))
			throw std::runtime_error("Bad offset order");

		auto const ref_begin(reference.cbegin());
		auto const output_start(ref_begin + output_start_pos);
		auto const output_end(ref_begin + output_end_pos);
		for (auto &kv : ref_haplotype_ptrs)
		{
			for (auto h_ptr : kv.second)
			{
				if (h_ptr)
				{
					auto &h(*h_ptr);
					if (h.current_pos != output_start_pos)
						throw std::runtime_error("Unexpected position");
					
					std::copy(output_start, output_end, std::ostream_iterator <char>(h.output_stream));
					h.current_pos = output_end_pos;
				}
			}
		}
	}
	
	
	void process_overlap_stack(
		overlap_stack_type &overlap_stack,
		size_t const var_pos,
		haplotype_ptr_map &ref_haplotype_ptrs,
		vector_type const &reference
	)
	{
		while (true)
		{
			// NOTE: for libc++, use p overlap_stack.c (not p overlap_stack).
			variant_overlap &vo(overlap_stack.top());
			if (var_pos < vo.end_pos)
				break;

			// The sequence was output up to vo's start_pos when it was added to the stack.
			// Check the current position and output from there.
			auto const output_start_pos(vo.current_pos);
			auto const output_end_pos(vo.end_pos);
			
			// Output reference from 5' direction up to vo.end_pos.
			output_reference(ref_haplotype_ptrs, reference, vo.current_pos, vo.end_pos);
			
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
							if (! (h.current_pos <= output_start_pos))
								throw std::runtime_error("Unexpected position");
							
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
				fill_streams(ref_haplotype_ptrs, fill_amt);
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
					auto &ref_ptrs(ref_haplotype_ptrs[sample_name]);
					
					for (size_t i(0), count(alt_ptrs.size()); i < count; ++i)
					{
						if (alt_ptrs[i] && ref_ptrs[i])
							throw std::runtime_error("Inconsistent haplotype pointers");
						
						// Use ADL.
						using std::swap;
						if (alt_ptrs[i])
							swap(alt_ptrs[i], ref_ptrs[i]);
					}
				}
			}

			// Update current_pos and heaviest path length with the result.
			if (1 < overlap_stack.size())
			{
				overlap_stack.pop();
				auto &previous_overlap(overlap_stack.top());
				previous_overlap.current_pos = output_end_pos;
				previous_overlap.heaviest_path_length += heaviest_path_length;
			}
			else
			{
				// If the handled variant_overlap was the last one, exit the loop.
				break;
			}
		}
	}


	void process_variants(
		vl::VariantCallFile &variant_file,
		variant_set const &skipped_variants,
		sample_set const &current_samples,
		haplotype_map &all_haplotypes,
		vector_type const &reference,
		std::string const &null_allele_seq
	)
	{
		haplotype_ptr_map ref_haplotype_ptrs;	// Haplotypes to which the reference sequence is to be output.
		overlap_stack_type overlap_stack;
		overlap_stack.emplace(0, 0, 0, 0);

		// All haplotypes initially have the reference sequence.
		for (auto &kv : all_haplotypes)
		{
			auto const &sample_name(kv.first);
			auto &haplotype_vector(kv.second);
			auto &haplotype_ptr_vector(ref_haplotype_ptrs[sample_name]);
			
			auto const count(haplotype_vector.size());
			haplotype_ptr_vector.resize(count);
			for (size_t i(0); i < count; ++i)
				haplotype_ptr_vector[i] = &haplotype_vector[i];
		}

		variant_file.reset();
		vl::Variant var(variant_file);
		size_t i(0);
		while (variant_file.getNextVariant(var))
		{
			if (0 != skipped_variants.count(var.id))
				continue;
			
			size_t const var_pos(zero_based_position_safe(var));

			// If var is beyond previous_variant.end_pos, handle the variants on the stack
			// until a containing variant is found or the bottom of the stack is reached.
			process_overlap_stack(overlap_stack, var_pos, ref_haplotype_ptrs, reference);
			
			// Use the previous variant's range to determine the output sequence.
			auto &previous_variant(overlap_stack.top());
			
			// Check that if var is before previous_variant.end_pos, it is also completely inside it.
			if (var_pos < previous_variant.end_pos && !(var_pos + var.ref.size() <= previous_variant.end_pos))
				throw std::runtime_error("Invalid variant inclusion");
			
			// Output reference from 5' direction up to var_pos.
			output_reference(ref_haplotype_ptrs, reference, previous_variant.current_pos, var_pos);
			
			// Add the amount output to the heaviest path length.
			previous_variant.heaviest_path_length += var_pos - previous_variant.current_pos;
			
			// Update current_pos to match the position up to which the sequence was output.
			// Also add the length to the heaviest path length.
			previous_variant.current_pos = var_pos;
			previous_variant.heaviest_path_length += var_pos - previous_variant.current_pos;
			
			// Find haplotypes that have the variant.
			alt_map alt_haplotypes;
			for (std::string const &sample : var.sampleNames)
			{
				auto &ref_ptrs(ref_haplotype_ptrs[sample]);
				
				auto const gt(var.getGenotype(sample));
				auto const decomposed(vl::decomposePhasedGenotype(gt));
				size_t chr_idx(0);
				for (auto const alt_idx : decomposed)
				{
					if (0 != alt_idx)
					{
						std::string const *alt(&null_allele_seq);
						if (vl::NULL_ALLELE != alt_idx)
							alt = &var.alt[alt_idx - 1];
						
						haplotype_ptr_map &alt_ptrs_by_sample(alt_haplotypes[*alt]);

						auto it(alt_ptrs_by_sample.find(sample));
						if (alt_ptrs_by_sample.end() == it)
						{
							it = alt_ptrs_by_sample.emplace(
								std::piecewise_construct,
								std::forward_as_tuple(sample),
								std::forward_as_tuple(ref_ptrs.size(), nullptr)
							).first;
						}
						auto &alt_ptrs(it->second);
						
						if ((alt_ptrs[chr_idx] && ref_ptrs[chr_idx]) || !(alt_ptrs[chr_idx] || ref_ptrs[chr_idx]))
							throw std::runtime_error("Inconsistent haplotype pointers");
						
						// Use ADL.
						using std::swap;
						swap(alt_ptrs[chr_idx], ref_ptrs[chr_idx]);
					}
					
					++chr_idx;
				}
			}
			
			// Create a new variant_overlap.
			variant_overlap overlap(var_pos, var_pos, var_pos + var.ref.size(), 0, alt_haplotypes);
			if (var_pos < previous_variant.end_pos)
			{
				// Add the current variant to the stack.
				overlap_stack.emplace(std::move(overlap));
			}
			else
			{
				// Replace the top of the stack with the current variant.
				// Use ADL.
				using std::swap;
				swap(overlap_stack.top(), overlap);
			}

			++i;
			if (0 == i % 100000)
				std::cerr << "Handled " << i << " variants…" << std::endl;
		}
		
		// Fill the remaining part with reference.
		process_overlap_stack(overlap_stack, reference.size(), ref_haplotype_ptrs, reference);
		
		{
			auto const ref_begin(reference.cbegin());
			auto const ref_end(reference.cend());
			
			for (auto &kv : all_haplotypes)
			{
				for (auto &h : kv.second)
					std::copy(ref_begin + h.current_pos, ref_end, std::ostream_iterator <char>(h.output_stream));
			}
		}
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	file_istream ref_fasta_stream;
	file_istream vcf_stream;

	open_file_for_reading(args_info.reference_arg, ref_fasta_stream);
	open_file_for_reading(args_info.variants_arg, vcf_stream);

	vl::VariantCallFile variant_file;
	vector_type reference;

	variant_file.open(vcf_stream);
	if (!variant_file.is_open())
	{
		std::cerr << "Opened the variant call file as a stream but was unable to create the file object." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	{
		// Read the reference file and place its contents into reference.
		typedef v2m::vector_source <vector_type> vector_source;

		struct callback {
			vector_type *reference;
			size_t i{0};

			callback(vector_type &ref):
				reference(&ref)
			{
			}

			void handle_sequence(
				std::string const &identifier,
				std::unique_ptr <vector_type> &seq,
				size_t const &seq_length, 
				vector_source &vs
			)
			{
				std::cerr << "Read sequence of length " << seq_length << std::endl;

				// Use ADL.
				using std::swap;
				seq->resize(seq_length);
				swap(*reference, *seq);
				vs.put_vector(seq);
			}

			void finish() {}
		};

		typedef v2m::fasta_reader <vector_source, callback> fasta_reader;

		vector_source vs(1, false);

		{
			// Approximate the size of the reference from the size of the FASTA file.
			int const fd(ref_fasta_stream->handle());
			struct stat sb;
			auto const st(fstat(fd, &sb));
			if (0 != st)
			{
				auto const err_str(strerror(errno));
				auto const msg(boost::str(boost::format("Unable to stat the reference file: %s") % err_str));
				throw(std::runtime_error(msg));
			}

			// Preallocate space for the reference.
			std::cerr << "Preallocating a vector of size " << sb.st_size << "…";
			std::unique_ptr <vector_type> vec_ptr;
			vs.get_vector(vec_ptr);
			vec_ptr->reserve(sb.st_size);
			vs.put_vector(vec_ptr);
			std::cerr << " done." << std::endl;
		}

		callback cb(reference);
		fasta_reader reader;

		std::cerr << "Reading reference FASTA into memory…" << std::endl;
		reader.read_from_stream(ref_fasta_stream, vs, cb);
	}

	{
		vl::Variant var(variant_file);

		if (!args_info.no_phasing_check_flag)
		{
			std::cerr << "Checking that the variant file is phased…" << std::endl;
			size_t i(0);
			variant_file.reset();
			while (variant_file.getNextVariant(var))
			{
				if (!var.isPhased())
				{
					std::cerr << "Variant " << var.id << " is not phased, cannot convert to FASTA." << std::endl;
					exit(EXIT_FAILURE);
				}

				++i;
				if (0 == i % 100000)
					std::cerr << "Handled " << i << " variants…" << std::endl;
			}
		}

		if (!args_info.no_unique_id_check_flag)
		{
			std::cerr << "Checking that the variant identifiers are unique…" << std::endl;
			size_t i(0);
			bool can_continue(true);
			std::set <std::string> identifiers;
			variant_file.reset();
			while (variant_file.getNextVariant(var))
			{
				auto const res(identifiers.insert(var.id));
				if (!res.second)
				{
					can_continue = false;
					std::cerr << "Found duplicate identifier '" << var.id << '\'' << std::endl;
				}

				++i;
				if (0 == i % 100000)
					std::cerr << "Handled " << i << " variants…" << std::endl;
			}
		}

		ploidy_map ploidy;
		if (!args_info.no_ploidy_check_flag)
			std::cerr << "Checking ploidy…" << std::endl;
		{
			size_t i(0);
			variant_file.reset();
			while (variant_file.getNextVariant(var))
			{
				ploidy_map current_ploidy;
				current_ploidy.clear();
				for (std::string const &sample : var.sampleNames)
				{
					auto const gt(var.getGenotype(sample));
					auto const decomposed(vl::decomposeGenotype(gt));
					for (auto const &kv : decomposed)
					{
						// first: allele, second: count
						current_ploidy[sample] += kv.second;
					}
				}
				
				if (! (
					ploidy.empty() ||
					std::equal(current_ploidy.cbegin(), current_ploidy.cend(), ploidy.cbegin(), ploidy.cend())
				))
				{
					std::cerr << "Ploidy changes in " << var.id << ", unable to continue." << std::endl;
					exit(EXIT_FAILURE);
				}
				
				// Use ADL.
				using std::swap;
				swap(ploidy, current_ploidy);

				if (args_info.no_ploidy_check_flag)
					break;

				++i;
				if (0 == i % 100000)
					std::cerr << "Handled " << i << " variants…" << std::endl;
			}
		}
		
		variant_set skipped_variants;
		{
			typedef boost::bimap <
				boost::bimaps::multiset_of <std::string>,
				boost::bimaps::multiset_of <std::string>
			> overlap_map;
			typedef boost::bimap <
				boost::bimaps::set_of <std::string>,
				boost::bimaps::list_of <size_t>
			> conflict_count_map;
			
			std::cerr << "Checking overlapping variants…" << std::endl;
			variant_file.reset();
			size_t last_position(0);
			std::map <size_t, std::string> end_positions;
			conflict_count_map conflict_counts;
			overlap_map bad_overlaps;
			size_t i(0);
			size_t conflict_count(0);
			while (variant_file.getNextVariant(var))
			{
				auto const pos(zero_based_position_safe(var));
				if (! (last_position <= pos))
					throw std::runtime_error("Positions not in increasing order");
				
				auto const end(pos + var.ref.size());
				auto it(end_positions.upper_bound(pos));
				auto const end_it(end_positions.cend());
				if (end_it == it)
				{
					auto const var_id(var.id);
					end_positions.emplace(end, var_id);
					goto loop_end;
				}
				
				// Found one or more possible bad overlaps. List all of them.
				do
				{
					// Check if the current variant is completely inside the possibly conflicting one.
					if (end <= it->first)
						break;
					
					++conflict_count;
					std::cerr << "Variant " << var.id << " conflicts with " << it->second << "." << std::endl;
					
					auto const res(bad_overlaps.insert(overlap_map::value_type(it->second, var.id)));
					if (false == res.second)
						throw std::runtime_error("Unable to insert");

					++conflict_counts.left[it->second];
					++conflict_counts.left[var.id];
					++it;
				} while (end_it != it);
				
				// Add the end position.
				{
					auto const var_id(var.id);
					end_positions.emplace(end, var_id);
				}

			loop_end:
				++i;
				if (0 == i % 100000)
					std::cerr << "Handled " << i << " variants…" << std::endl;
			}
			
			
			// Remove conflicting variants starting from the one with the highest score.
			conflict_counts.right.sort();
			for (auto const &kv : boost::adaptors::reverse(conflict_counts.right))
			{
				auto const count(kv.first);
				auto const &var_id(kv.second);
				
				skipped_variants.insert(var_id);
				bad_overlaps.left.erase(var_id);
				bad_overlaps.right.erase(var_id);
				
				if (bad_overlaps.size() == 0)
					break;
			}
			
			if (bad_overlaps.size() != 0)
				throw std::runtime_error("Unable to remove all conflicting variants");
			
			auto const skipped_count(skipped_variants.size());
			if (0 == skipped_count)
				std::cerr << "Found no conflicting variants." << std::endl;
			else
			{
				std::cerr << "Found " << conflict_count << " conflicts in total." << std::endl;
				std::cerr << "Skipping the following conflicting variants:" << std::endl;
				for (auto const &id : skipped_variants)
					std::cerr << '\t' << id << std::endl;

				std::cerr << "Number of variants to be skipped: " << skipped_variants.size() << std::endl;
			}
		}

		std::cerr << "Generating haplotype sequences…" << std::endl;
		auto const &all_samples(variant_file.sampleNames);
		auto const sample_count(all_samples.size());
		size_t const chunk_size(200);
		size_t const total_rounds(std::ceil(1.0 * sample_count / chunk_size));
		sample_set current_samples;
		std::string const null_allele_seq("N"); // FIXME: use an argument.

		haplotype_map haplotypes;
		for (size_t i(0); i < sample_count; i += chunk_size)
		{
			std::cerr << "Round " << (i + 1) << '/' << total_rounds << std::endl;

			// Update current_samples.
			current_samples.clear();
			auto begin(all_samples.cbegin() + i);
			auto const end(all_samples.cbegin() + std::min(i + chunk_size, sample_count));
			std::copy(begin, end, std::inserter(current_samples, current_samples.end()));

			open_files_for_writing(current_samples, haplotypes, ploidy, args_info.overwrite_flag);

			// Move to the beginning of the variant file.
			process_variants(variant_file, skipped_variants, current_samples, haplotypes, reference, null_allele_seq);

		}
	}

	return EXIT_SUCCESS;
}
