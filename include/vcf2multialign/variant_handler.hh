/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_HANDLER_HH
#define VCF2MULTIALIGN_VARIANT_HANDLER_HH

#include <dispatch/dispatch.h>
#include <libbio/assert.hh>
#include <libbio/dispatch_fn.hh>
#include <libbio/file_handling.hh>
#include <libbio/vcf_reader.hh>
#include <map>
#include <stack>
#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/variant_buffer.hh>


namespace vcf2multialign {

	enum { REF_SAMPLE_NUMBER = 0 };
	
	struct haplotype
	{
		size_t current_pos{0};
	};
	
	template <typename t_haplotype>
	using haplotype_map = std::map <
		std::size_t,				// Sample (column) number
		std::vector <t_haplotype>	// All haplotype sequences
	>;
	
	template <typename t_haplotype>
	using haplotype_ptr_map = std::map <
		std::size_t,				// Sample (column) number
		std::vector <t_haplotype *>	// Haplotype sequences by chromosome index
	>;
	
	template <typename t_haplotype>
	using alt_map = std::map <
		std::string,				// ALT
		haplotype_ptr_map <t_haplotype>
	>;
	
	
	struct skipped_sample
	{
		std::size_t	sample_no{0};
		uint8_t		alt_idx{0};
		uint8_t		chr_idx{0};
		
		skipped_sample(std::size_t const sample_no_, uint8_t const alt_idx_, uint8_t const chr_idx_):
			sample_no(sample_no_),
			alt_idx(alt_idx_),
			chr_idx(chr_idx_)
		{
		}
	};
	
	
	template <typename t_haplotype>
	struct variant_overlap
	{
		typedef alt_map <t_haplotype> alt_map_type;
		
		size_t			start_pos{0};
		size_t			current_pos{0};
		size_t			end_pos{0};
		size_t			heaviest_path_length{0};
		size_t			lineno{0};
		alt_map_type	alt_haplotypes;
		
		variant_overlap(
			size_t const start_pos_,
			size_t const current_pos_,
			size_t const end_pos_,
			size_t const heaviest_path_length_,
			size_t const lineno_
		):
			start_pos(start_pos_),
			current_pos(current_pos_),
			end_pos(end_pos_),
			heaviest_path_length(heaviest_path_length_),
			lineno(lineno_)
		{
			libbio::always_assert(start_pos <= end_pos, "Bad offset order");
		}
		
		variant_overlap(
			size_t const start_pos_,
			size_t const current_pos_,
			size_t const end_pos_,
			size_t const heaviest_path_length_,
			size_t const lineno_,
			alt_map_type &alts
		):
			start_pos(start_pos_),
			current_pos(current_pos_),
			end_pos(end_pos_),
			heaviest_path_length(heaviest_path_length_),
			lineno(lineno_),
			alt_haplotypes(std::move(alts))
		{
			libbio::always_assert(start_pos <= end_pos, "Bad offset order");
		}
	};
	
	
	class variant_handler_base
	{
	protected:
		std::set <size_t>								m_valid_alts;
		error_logger									*m_error_logger{};
		sv_handling										m_sv_handling_method{};
		
	protected:
		variant_handler_base() = default;
		
		variant_handler_base(
			sv_handling const sv_handling_method,
			error_logger &error_logger
		):
			m_error_logger(&error_logger),
			m_sv_handling_method(sv_handling_method)
		{
		}
		
		bool check_alt_seq(std::string const &alt) const;
		void fill_valid_alts(libbio::variant const &var);
	};
	
	
	template <typename t_haplotype, typename t_delegate>
	class variant_handler : public variant_handler_base, public variant_buffer_delegate
	{
		friend t_delegate;
		
	public:
		typedef haplotype_map <t_haplotype>				haplotype_map_type;
		
	protected:
		typedef variant_overlap <t_haplotype>			variant_overlap_type;
		typedef haplotype_ptr_map <t_haplotype>			haplotype_ptr_map_type;
		typedef alt_map <t_haplotype>					alt_map_type;
		typedef std::stack <variant_overlap_type>		overlap_stack_type;
		typedef std::vector <size_t>					sample_number_vector;

	protected:
		libbio::dispatch_ptr <dispatch_queue_t>			m_main_queue{};
		libbio::dispatch_ptr <dispatch_queue_t>			m_parsing_queue{};
		t_delegate										*m_delegate{};
		
		vector_type	const								*m_reference{};
		
		variant_buffer									m_variant_buffer;
		overlap_stack_type								m_overlap_stack;
		
		variant_set const								*m_skipped_variants{};
		variant_set										m_overlapping_alts{};
		haplotype_ptr_map_type							m_ref_haplotype_ptrs;			// Haplotypes to which the reference sequence is to be output.
		std::vector <skipped_sample>					m_skipped_samples;				// In current variant.
		std::map <uint8_t, sample_count>				m_counts_by_alt;				// In current variant.
		sample_count									m_non_ref_totals;				// In current variant.
		
		haplotype_map_type								*m_all_haplotypes{};
		alt_map_type									m_alt_haplotypes;
		
		std::string const								*m_null_allele_seq{};
		std::size_t										m_i{0};
		
	public:
		variant_handler(
			libbio::dispatch_ptr <dispatch_queue_t> const &main_queue,
			libbio::dispatch_ptr <dispatch_queue_t> const &parsing_queue,
			libbio::vcf_reader &vcf_reader_,
			vector_type const &reference,
			sv_handling const sv_handling_method,
			variant_set const &skipped_variants,
			std::string const &null_allele,
			error_logger &error_logger,
			t_delegate &delegate
		):
			variant_handler_base(sv_handling_method, error_logger),
			m_main_queue(main_queue),
			m_parsing_queue(parsing_queue),
			m_delegate(&delegate),
			m_reference(&reference),
			m_variant_buffer(vcf_reader_, main_queue, *this),
			m_skipped_variants(&skipped_variants),
			m_null_allele_seq(&null_allele)
		{
		}
		
		variant_handler() = default;
		
	public:
		void process_variants(haplotype_map_type &haplotypes);
		variant_buffer &get_variant_buffer() { return m_variant_buffer; }

	protected:
		virtual void handle_variant(libbio::variant &var);
		virtual void finish();

		void fill_streams(haplotype_ptr_map_type &haplotypes, size_t const fill_amt) const;
		void output_reference(std::size_t const output_start_pos, std::size_t const output_end_pos);
		std::size_t process_overlap_stack(size_t const var_pos);
	};
	
	
	// Fill the streams with '-'.
	template <typename t_haplotype, typename t_delegate>
	void variant_handler <t_haplotype, t_delegate>::fill_streams(haplotype_ptr_map_type &haplotypes, size_t const fill_amt) const
	{
		for (auto &kv : haplotypes)
		{
			for (auto h_ptr : kv.second)
			{
				if (h_ptr)
					h_ptr->fill_with_dashes(fill_amt);
			}
		}
	}
	
	
	// Fill the streams with reference.
	template <typename t_haplotype, typename t_delegate>
	void variant_handler <t_haplotype, t_delegate>::output_reference(std::size_t const output_start_pos, std::size_t const output_end_pos)
	{
		namespace lb = libbio;
		
		if (output_start_pos == output_end_pos)
			return;
		
		lb::always_assert(output_start_pos < output_end_pos, "Bad offset order");
		
		char const *ref_begin(m_reference->data());
		auto const output_start(ref_begin + output_start_pos);
		for (auto &kv : m_ref_haplotype_ptrs)
		{
			for (auto h_ptr : kv.second)
			{
				if (h_ptr)
				{
					auto &h(*h_ptr);
					lb::always_assert(h.current_pos == output_start_pos, "Unexpected position");
					
					h.write(output_start, output_end_pos - output_start_pos);
					h.current_pos = output_end_pos;
				}
			}
		}
	}
	
	
	template <typename t_haplotype, typename t_delegate>
	std::size_t variant_handler <t_haplotype, t_delegate>::process_overlap_stack(size_t const var_pos)
	{
		namespace lb = libbio;
		
		std::size_t retval(0);
		while (true)
		{
			// NOTE: for libc++, use p overlap_stack.c in the debugger (not p overlap_stack).
			auto &vo(m_overlap_stack.top());
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
							lb::always_assert(h.current_pos <= output_start_pos, "Unexpected position");
							
							h.append(alt);
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
						lb::always_assert(! (alt_ptrs[i] && ref_ptrs[i]), "Inconsistent haplotype pointers");
						
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
	
	
	template <typename t_haplotype, typename t_delegate>
	void variant_handler <t_haplotype, t_delegate>::handle_variant(libbio::variant &var)
	{
		namespace lb = libbio;
		
		m_skipped_samples.clear();
		m_counts_by_alt.clear();
		m_non_ref_totals.reset();
		
		auto const lineno(var.lineno());
		if (0 != m_skipped_variants->count(lineno))
			return;
		
		// Preprocess the ALT field to check that it can be handled.
		fill_valid_alts(var);
		if (m_valid_alts.empty())
			return;
		
		size_t const var_pos(var.zero_based_pos());
		lb::always_assert(var_pos < m_reference->size(), [this, lineno](){
			std::cerr
			<< "Variant position on line " << lineno
			<< " greater than reference length (" << m_reference->size() << ")."
			<< std::endl;
		});
		
		auto const var_ref(var.ref());
		auto const var_ref_size(var_ref.size());
		auto const var_alts(var.alts());
		auto const var_alt_sv_types(var.alt_sv_types());
		
		// If var is beyond previous_variant.end_pos, handle the variants on the stack
		// until a containing variant is found or the bottom of the stack is reached.
		process_overlap_stack(var_pos);
		m_delegate->variant_handler_did_process_overlap_stack(*this);
		
		// Use the previous variant's range to determine the output sequence.
		auto &previous_variant(m_overlap_stack.top());
		
		// Check that if var is before previous_variant.end_pos, it is also completely inside it.
		lb::always_assert(
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
		std::string const empty_alt("");
		for (auto const alt_idx : m_valid_alts)
		{
			switch (var_alt_sv_types[alt_idx - 1])
			{
				case lb::sv_type::NONE:
				{
					auto const &alt_str(var_alts[alt_idx - 1]);
					m_alt_haplotypes[alt_str];
					break;
				}
				
				case lb::sv_type::DEL:
				case lb::sv_type::DEL_ME:
					m_alt_haplotypes[empty_alt];
					break;
				
				default:
					lb::fail("Unexpected structural variant type.");
					break;
			}
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
				lb::always_assert(0 == chr_idx || is_phased, "Variant file not phased");
				
				if (0 != alt_idx && 0 != m_valid_alts.count(alt_idx))
				{
					auto &ref_ptrs(m_ref_haplotype_ptrs.find(sample_no)->second); // Has nodes for every sample_no.
					
					std::string const *alt_ptr{m_null_allele_seq};
					if (lb::NULL_ALLELE != alt_idx)
					{
						switch (var_alt_sv_types[alt_idx - 1])
						{
							case lb::sv_type::NONE:
								alt_ptr = &var_alts[alt_idx - 1];
								break;
								
							case lb::sv_type::DEL:
							case lb::sv_type::DEL_ME:
								alt_ptr = &empty_alt;
								break;
								
							default:
								lb::fail("Unexpected structural variant type.");
								break;
						}
					}
					
					auto &alt_ptrs_by_sample(m_alt_haplotypes[*alt_ptr]);
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
		variant_overlap_type overlap(var_pos, var_pos, var_end, 0, lineno, m_alt_haplotypes);
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
	
	
	template <typename t_haplotype, typename t_delegate>
	void variant_handler <t_haplotype, t_delegate>::finish()
	{
		// Fill the remaining part with reference.
		m_error_logger->flush();
		std::cerr << "Filling with the reference…" << std::endl;
		auto const ref_size(m_reference->size());
		auto const output_end_pos(process_overlap_stack(ref_size));
		m_delegate->variant_handler_did_process_overlap_stack(*this);
		
		char const *ref_begin(m_reference->data());
		for (auto &kv : *m_all_haplotypes)
		{
			for (auto &h : kv.second)
			{
				auto const output_len(ref_size - h.current_pos);
				h.write(ref_begin + h.current_pos, output_len);
			}
		}
		
		m_delegate->variant_handler_did_finish(*this);
	}
	
	
	template <typename t_haplotype, typename t_delegate>
	void variant_handler <t_haplotype, t_delegate>::process_variants(haplotype_map_type &all_haplotypes)
	{
		namespace lb = libbio;
		
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
		reader.set_parsed_fields(lb::vcf_field::ALL);
		
		lb::dispatch_caller caller(&m_variant_buffer);
		caller.template async <&variant_buffer::read_input>(*m_parsing_queue);
	}
}

#endif
