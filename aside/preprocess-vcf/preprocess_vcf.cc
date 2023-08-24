/*
 * Copyright (c) 2019–2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <libbio/progress_indicator.hh>
#include <libbio/utility.hh>
#include <vcf2multialign/preprocess/preprocess_logger.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/utility/dispatch_cli_runner.hh>
#include <vcf2multialign/utility/dispatch_exit_guard.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "preprocess_vcf.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	class progress_indicator_delegate final : public lb::progress_indicator_delegate
	{
	protected:
		v2m::variant_partitioner const * const	m_partitioner{};
		
	public:
		progress_indicator_delegate(v2m::variant_partitioner const &partitioner):
			m_partitioner(&partitioner)
		{
		}
		
		virtual std::size_t progress_step_max() const override { return 0; }
		virtual std::size_t progress_current_step() const override { return m_partitioner->processed_count(); }
		virtual void progress_log_extra() const override {};
	};
}


namespace vcf2multialign {
	
	struct ploidy_result
	{
		std::size_t		donor_count{};
		std::size_t 	chromosome_copy_count{};
		
		bool is_valid() const { return 0 < donor_count; }
	};
	
	
	class cut_position_processor final :	public dispatch_cli_runner,
											public vcf_processor,
											public output_stream_controller,
											public reference_controller,
											public preprocess_logger
	{
	protected:
		std::vector <std::string> const	m_field_names_for_filter_if_set;
		std::string	const				m_chr_name;
		std::size_t const				m_minimum_subgraph_distance{};
		
	public:
		cut_position_processor(
			std::vector <std::string> &&field_names_for_filter_if_set,
			std::string const &chr_name,
			std::size_t const minimum_subgraph_distance
		):
			dispatch_cli_runner(),
			vcf_processor(),
			output_stream_controller(),
			m_field_names_for_filter_if_set(std::move(field_names_for_filter_if_set)),
			m_chr_name(chr_name),
			m_minimum_subgraph_distance(minimum_subgraph_distance)
		{
		}
		
		ploidy_result check_donor_and_chromosome_copy_count();
		
	protected:
		void do_work() override;
	};
	
	
	void cut_position_processor::do_work()
	{
		// Set up the processor.
		variant_partitioner partitioner(
			*this,
			this->m_vcf_reader,
			this->m_reference,
			this->m_chr_name,
			this->m_minimum_subgraph_distance
		);
		
		preprocessing_result result;
		
		// Partition.
		progress_indicator_delegate indicator_delegate(partitioner);
		this->install_progress_indicator();
		
		this->progress_indicator().log_with_counter(lb::copy_time() + "Processing the variants…", indicator_delegate);
		partitioner.partition(m_field_names_for_filter_if_set, result, false);
		this->end_logging();
		this->uninstall_progress_indicator();
		
		dispatch_async(dispatch_get_main_queue(), ^{
			lb::log_time(std::cerr);
			std::cerr
				<< "Done. Handled variants: " << result.handled_line_numbers.size()
				<< " Maximum segment size: " << result.max_segment_size
				<< " Cut position count: " << result.positions.size()
				<< " Suitable variants: " << this->variants_passing_checks()
				<< " Chromosome ID mismatches: " << this->chrom_id_mismatches()
				<< '\n';
		});
		
		// Output.
		cereal::PortableBinaryOutputArchive archive(this->m_output_stream);
		archive(result);
		this->m_output_stream << std::flush;
	}
	
	
	void preprocess_vcf(
		char const *reference_path,
		char const *variant_file_path,
		char const *output_path,
		char const *log_path,
		char const *reference_seq_name,
		char const *chr_name,
		std::vector <std::string> &&field_names_for_filter_if_set,
		std::size_t const minimum_subgraph_distance,
		bool const should_overwrite_files
	)
	{
		// main() calls dispatch_main() if this function does not throw or call std::exit or something similar.
		try
		{
			// Since the processor has Boost’s streams, it cannot be moved. Hence the use of a pointer.
			typedef v2m::dispatch_exit_guard_helper <cut_position_processor> wrapped_processor_type;
			
			auto processor_ptr(std::make_unique <wrapped_processor_type>(std::move(field_names_for_filter_if_set), chr_name, minimum_subgraph_distance));
			auto &processor(processor_ptr->value);
			
			// These will eventually call std::exit if the file in question cannot be opened.
			processor.open_variants_file(variant_file_path);
			processor.read_reference(reference_path, reference_seq_name);
			processor.open_output_file(output_path, should_overwrite_files);
			if (log_path)
				processor.open_log_file(log_path, should_overwrite_files);
			
			processor.prepare_reader();
			
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
