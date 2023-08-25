/*
 * Copyright (c) 2019–2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <libbio/progress_indicator.hh>
#include <tuple>
#include <unistd.h>
#include <vcf2multialign/preprocess/preprocess_logger.hh>
#include <vcf2multialign/utility/can_handle_variant_alts.hh>
#include <vcf2multialign/utility/dispatch_cli_runner.hh>
#include <vcf2multialign/utility/dispatch_exit_guard.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/variant_format.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace {
	
	class progress_indicator_delegate final : public lb::progress_indicator_delegate
	{
	protected:
		std::atomic_size_t	m_counter{};
		
	public:
		progress_indicator_delegate() = default;
		void increment() { m_counter.fetch_add(1, std::memory_order_relaxed); }
		virtual std::size_t progress_step_max() const override { return 0; }
		virtual std::size_t progress_current_step() const override { return m_counter.load(std::memory_order_relaxed); }
		virtual void progress_log_extra() const override {};
	};
}


namespace vcf2multialign {

	struct variant_statistics
	{
		std::size_t	included_variant{};
		std::size_t	ref_mismatch{};
		std::size_t	filtered_variant{};
		std::size_t	no_suitable_alts{};
		std::size_t	gt_equals_0{};
		std::size_t	overlapping_variant{};
		std::size_t	unhandled_structural_variant{};

		variant_statistics() = default;
	};

	std::ostream &operator<<(std::ostream &os, variant_statistics const &stats)
	{
		os << "Included variants:             " << stats.included_variant             << '\n';
		os << "Reference mismatch:            " << stats.ref_mismatch                 << '\n';
		os << "Filtered variants:             " << stats.filtered_variant             << '\n';
		os << "Overlapping variants:          " << stats.overlapping_variant          << '\n';
		os << "No suitable ALTs:              " << stats.no_suitable_alts             << '\n';
		os << "GT equals 0:                   " << stats.gt_equals_0                  << '\n';
		os << "Unhandled structural variants: " << stats.unhandled_structural_variant << '\n';
		return os;
	}

	
	class unaligned_processor final :	public dispatch_cli_runner,
										public vcf_processor,
										public output_stream_controller,
										public reference_controller,
										public preprocess_logger
	{
	protected:
		std::vector <std::string> const	m_excluded_filters;
		std::vector <std::string> const	m_field_names_for_filter_if_set;
		std::string const				m_chromosome_name;
		std::size_t const				m_sample_idx;
		std::uint16_t const				m_chromosome_copy_idx;
		
		
	public:
		unaligned_processor(
			std::vector <std::string> &&excluded_filters,
			std::vector <std::string> &&field_names_for_filter_if_set,
			std::string const &chromosome_name,
			std::size_t const sample_idx,
			std::uint16_t const chromosome_copy_idx
		):
			dispatch_cli_runner(),
			vcf_processor(),
			output_stream_controller(),
			reference_controller(),
			preprocess_logger(),
			m_excluded_filters(std::move(excluded_filters)),
			m_field_names_for_filter_if_set(std::move(field_names_for_filter_if_set)),
			m_chromosome_name(chromosome_name),
			m_sample_idx(sample_idx),
			m_chromosome_copy_idx(chromosome_copy_idx)
		{
		}
		
		std::vector <std::string> const &excluded_filters() const { return m_excluded_filters; }
		std::vector <std::string> const &field_names_for_filter_if_set() const { return m_field_names_for_filter_if_set; }
		
	protected:
		void do_work() override;
		void output_variants(variant_statistics &statistics, progress_indicator_delegate &progress_delegate);
	};
	
	
	void unaligned_processor::do_work()
	{
		progress_indicator_delegate indicator_delegate;
		
		this->progress_indicator().log_with_counter(lb::copy_time() + "Processing the variants…", indicator_delegate);
		
		variant_statistics statistics;
		output_variants(statistics, indicator_delegate);
		
		this->end_logging();
		this->uninstall_progress_indicator();
		
		dispatch_async(dispatch_get_main_queue(), ^{ 
			lb::log_time(std::cerr);
			std::cerr << "Done.\n";
			std::cerr << statistics;
		});
		
		this->m_output_stream << std::flush;
	}
	
	
	void unaligned_processor::output_variants(variant_statistics &statistics, progress_indicator_delegate &indicator_delegate)
	{
		auto &reader(this->vcf_reader());
		reader.set_parsed_fields(vcf::field::ALL);
		
		// Get the field descriptors needed for accessing the values.
		auto const *end_field(reader.get_end_field_ptr());
		
		// Determine the fields used for filtering.
		std::vector <vcf::info_field_base *> filter_by_assigned;
		{
			auto const &fields(reader.info_fields());
			for (auto const &name : m_field_names_for_filter_if_set)
			{
				auto const it(fields.find(name));
				if (fields.end() == it)
				{
					this->variant_processor_no_field_for_identifier(name);
					continue;
				}
				
				filter_by_assigned.emplace_back(it->second.get());
			}
		}
		
		bool should_continue(false);
		std::size_t prev_end_pos{};
		std::size_t overlap_end{};
		reader.parse(
			[
				this,
				end_field,
				&prev_end_pos,
				&filter_by_assigned,
				&statistics,
				&indicator_delegate
			](vcf::transient_variant const &var) -> bool
			{
				auto const lineno(var.lineno());
				auto const var_pos(var.zero_based_pos());
				
				// Check the chromosome name.
				if (var.chrom_id() != m_chromosome_name)
					goto end;
				
				if (! (var_pos < m_reference.size()))
				{
					this->variant_processor_found_variant_with_position_greater_than_reference_length(var);
					return false;
				}
				
				if (!can_handle_variant_alts(var))
				{
					this->variant_processor_found_variant_with_no_suitable_alts(var);
					++statistics.no_suitable_alts;
					goto end;
				}
				
				// Filter.
				{
					auto const &filters(var.filters());
					for (auto const *filter : filters)
					{
						auto const &filter_id(filter->get_id());
						if (std::binary_search(m_excluded_filters.begin(), m_excluded_filters.end(), filter_id))
						{
							++statistics.filtered_variant;
							goto end;
						}
					}
				}

				for (auto const *field_ptr : filter_by_assigned)
				{
					if (field_ptr->has_value(var))
					{
						this->variant_processor_found_filtered_variant(var, *field_ptr);
						++statistics.filtered_variant;
						goto end;
					}
				}
				
				// Compare the REF column against the reference sequence.
				{
					auto const &ref_col(var.ref());
					std::string_view const ref_sub(m_reference.data() + var_pos, ref_col.size());
					if (ref_col != ref_sub)
					{
						this->variant_processor_found_variant_with_ref_mismatch(var, ref_sub);
						++statistics.ref_mismatch;
						goto end;
					}
				}
				
				// Variant passes the checks, handle it.
				{
					libbio_always_assert_lt(m_sample_idx, var.samples().size());
					auto const &sample(var.samples()[m_sample_idx]);
					auto const *gt_field(get_variant_format(var).gt);
					auto const &gt((*gt_field)(sample));
					libbio_always_assert_lt(m_chromosome_copy_idx, gt.size());
					auto const alt_idx(gt[m_chromosome_copy_idx].alt);
					if (0 == alt_idx)
					{
						++statistics.gt_equals_0;
						goto end;
					}
					else
					{
						if (! (prev_end_pos <= var_pos))
						{
							// FIXME: log the overlap.
							//this->variant_processor_found_overlapping_variant(var, prev_end_pos);
							++statistics.overlapping_variant;
							goto end;
						}
						
						auto const &alt(var.alts()[alt_idx - 1]);
						switch (alt.alt_sv_type)
						{
							// Handled ALTs:
							case vcf::sv_type::NONE:
							case vcf::sv_type::UNKNOWN:
							case vcf::sv_type::DEL:
							case vcf::sv_type::DEL_ME:
								break;
							default:
								// FIXME: log the ALT.
								//this->variant_processor_unhandled_alt(var, alt_idx);
								++statistics.unhandled_structural_variant;
								goto end;
						}
						
						std::string_view const ref_sub(m_reference.data() + prev_end_pos, var_pos - prev_end_pos);
						this->m_output_stream << ref_sub;
						
						if (vcf::sv_type::NONE == alt.alt_sv_type)
						{
							auto const &alt_seq(alt.alt);
							this->m_output_stream << alt_seq;
						}
						
						auto const var_end(vcf::variant_end_pos(var, *end_field));
						prev_end_pos = var_end;
						++statistics.included_variant;
					}

					goto end;
				}
				
			end:
				indicator_delegate.increment();
				return true;
			}
		);
		
		std::string_view const ref_sub(m_reference.data() + prev_end_pos, m_reference.size() - prev_end_pos);
		this->m_output_stream << ref_sub;
	}
	
	
	void vcf_to_unaligned(
		char const *reference_path,
		char const *variants_file_path,
		char const *output_file_path,
		char const *log_path,
		char const *reference_seq_name,
		char const *chr_name,
		std::vector <std::string> &&excluded_filters,
		std::vector <std::string> &&field_names_for_filter_if_set,
		std::size_t sample_idx,
		std::uint16_t chr_copy_idx,
		bool const should_overwrite_files
	)
	{
		// main() calls dispatch_main() if this function does not throw or call std::exit or something similar.
		try
		{
			std::sort(excluded_filters.begin(), excluded_filters.end());

			// Since the processor has Boost’s streams, it cannot be moved. Hence the use of a pointer.
			typedef dispatch_exit_guard_helper <unaligned_processor> wrapped_processor_type;
			auto processor_ptr(std::make_unique <wrapped_processor_type>(
				std::move(excluded_filters),
				std::move(field_names_for_filter_if_set),
				chr_name,
				sample_idx,
				chr_copy_idx
			));
			auto &processor(processor_ptr->value);
			
			// These will eventually call std::exit if the file in question cannot be opened.
			processor.open_variants_file(variants_file_path);
			processor.read_reference(reference_path, reference_seq_name);
			processor.open_output_file(output_file_path, should_overwrite_files);
			if (log_path)
				processor.open_log_file(log_path, should_overwrite_files);
			
			processor.prepare_reader();

			// Check whether the given field IDs are listed in VCF headers.
			{
				auto const &reader(processor.vcf_reader());
				auto const &metadata(reader.metadata());
				auto const &filters(metadata.filter());
				auto const &info_fields(metadata.info());
				auto const &excluded_filters(processor.excluded_filters());
				auto const &field_names_for_filter_if_set(processor.field_names_for_filter_if_set());

				for (auto const &filter_id : excluded_filters)
				{
					if (filters.end() == filters.find(filter_id))
						std::cerr << "WARNING: Filter with ID ‘" << filter_id << "’ not found in VCF headers.\n";
				}

				for (auto const &info_field_id : field_names_for_filter_if_set)
				{
					if (info_fields.end() == info_fields.find(info_field_id))
						std::cerr << "WARNING: Info field with ID ‘" << info_field_id << "’ not found in VCF headers.\n";
				}
			}
			
			// Run in background in order to be able to update a progress bar.
			lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
				[
					processor_ptr = std::move(processor_ptr)
				](){
					processor_ptr->value.run(true);
				}
			);
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			log_assertion_failure_exception(exc);
			throw exc;
		}
	}
}


int main(int argc, char **argv)
{
#ifndef NDEBUG
	std::cerr << "Assertions have been enabled." << std::endl;
#endif
	
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	if (! (0 <= args_info.sample_index_arg))
	{
		std::cerr << "Sample index must be greater than or equal to zero.\n";
		std::exit(EXIT_FAILURE);
	}
	
	if (! (0 <= args_info.chromosome_copy_index_arg))
	{
		std::cerr << "Chromosome copy index must be greater than or equal to zero.\n";
		std::exit(EXIT_FAILURE);
	}

	std::vector <std::string> excluded_filters(args_info.exclude_filter_given);
	for (std::size_t i(0); i < args_info.exclude_filter_given; ++i)
		excluded_filters[i] = args_info.exclude_filter_arg[i];
	
	std::vector <std::string> field_names_for_filter_if_set(args_info.filter_field_set_given);
	for (std::size_t i(0); i < args_info.filter_field_set_given; ++i)
		field_names_for_filter_if_set[i] = args_info.filter_field_set_arg[i];
	
	v2m::vcf_to_unaligned(
		args_info.reference_arg,
		args_info.variants_arg,
		args_info.output_arg,
		args_info.log_arg,
		args_info.reference_sequence_given ? args_info.reference_sequence_arg : nullptr,
		args_info.chromosome_given ? args_info.chromosome_arg : nullptr,
		std::move(excluded_filters),
		std::move(field_names_for_filter_if_set),
		args_info.sample_index_arg,
		args_info.chromosome_copy_index_arg,
		args_info.overwrite_flag
	);
	
	dispatch_main();
	// Not reached.
	return EXIT_SUCCESS;
}