/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <libbio/progress_indicator.hh>
#include <libbio/utility.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/utility/check_ploidy.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/utility/progress_indicator_manager.hh>
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
	
	
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
}


namespace vcf2multialign {
	
	class cut_position_processor :	public vcf_processor,
									public output_stream_controller,
									public reference_controller,
									public progress_indicator_manager
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
			vcf_processor(),
			output_stream_controller(),
			progress_indicator_manager(),
			m_field_names_for_filter_if_set(std::move(field_names_for_filter_if_set)),
			m_chr_name(chr_name),
			m_minimum_subgraph_distance(minimum_subgraph_distance)
		{
		}
		
		void process_and_output(std::size_t const donor_count, std::size_t const chr_count);
	};
	
	
	void cut_position_processor::process_and_output(std::size_t const donor_count, std::size_t const chr_count)
	{
		try
		{
			// Set up the processor.
			logging_variant_processor_delegate processor_delegate;
			variant_partitioner partitioner(
				processor_delegate,
				this->m_vcf_reader,
				this->m_reference,
				this->m_chr_name,
				donor_count,
				chr_count,
				this->m_minimum_subgraph_distance
			);
			
			cut_position_list cut_positions;
			cut_positions.donor_count = donor_count;
			cut_positions.chr_count = chr_count;
			
			// Partition.
			progress_indicator_delegate indicator_delegate(partitioner);
			this->install_progress_indicator();
			
			this->progress_indicator().log_with_counter(lb::copy_time() + "Processing the variants…", indicator_delegate);
			partitioner.partition(m_field_names_for_filter_if_set, cut_positions);
			this->end_logging();
			this->uninstall_progress_indicator();
			
			dispatch_async_main(^{
				lb::log_time(std::cerr);
				std::cerr << "Done. Maximum segment size: " << cut_positions.max_segment_size << " Cut position count: " << cut_positions.positions.size() << '\n';
			});
			
			// Output.
			cereal::PortableBinaryOutputArchive archive(this->m_output_stream);
			archive(cut_positions);
			this->m_output_stream << std::flush;
			
			this->finish();
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			this->log_assertion_failure_and_exit(exc);
		}
		catch (std::exception const &exc)
		{
			this->log_exception_and_exit(exc);
		}
		catch (...)
		{
			this->log_unknown_exception_and_exit();
		}
	}
	
	
	void preprocess_vcf(
		char const *reference_path,
		char const *variant_file_path,
		char const *output_path,
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
			auto processor_ptr(std::make_unique <cut_position_processor>(std::move(field_names_for_filter_if_set), chr_name, minimum_subgraph_distance));
			auto &processor(*processor_ptr);
			
			// These will eventually call std::exit if the file in question cannot be opened.
			processor.open_variants_file(variant_file_path);
			processor.read_reference(reference_path, reference_seq_name);
			processor.open_output_file(output_path, should_overwrite_files);
			
			processor.prepare_reader();
			
			// Check the ploidy.
			ploidy_map ploidy;
			check_ploidy(processor.vcf_reader(), ploidy);
			auto const donor_count(ploidy.size());
			if (!donor_count)
			{
				std::cerr << "WARNING: No donors found." << std::endl;
				std::exit(EXIT_SUCCESS);
			}
			auto const chr_count(ploidy.begin()->second);
			
			// Run in background in order to be able to update a progress bar.
			lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
				[
					processor_ptr = std::move(processor_ptr),
				 	chr_count,
				 	donor_count
				](){
					processor_ptr->process_and_output(donor_count, chr_count);
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
