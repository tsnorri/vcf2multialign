/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <libbio/progress_indicator.hh>
#include <vcf2multialign/graph/variant_graph_generator.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/utility/progress_indicator_manager.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "create_variant_graph.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	class variant_graph_processor; // Fwd.
}


namespace {
	
	class progress_indicator_delegate final : public lb::progress_indicator_delegate
	{
	protected:
		v2m::variant_graph_processor const * const	m_processor{};
		
	public:
		progress_indicator_delegate(v2m::variant_graph_processor const &processor):
			m_processor(&processor)
		{
		}
		
		virtual void progress_log_extra() const override {};
		inline virtual std::size_t progress_step_max() const override;
		inline virtual std::size_t progress_current_step() const override;
	};
}


namespace vcf2multialign {
	
	class variant_graph_processor :	public vcf_processor,
									public output_stream_controller,
									public reference_controller,
									public progress_indicator_manager,
									public variant_graph_generator_delegate
	{
		friend progress_indicator_delegate;
		
	protected:
		preprocessing_result	m_preprocessing_result;
		variant_graph_generator	m_generator;
		
	public:
		variant_graph_processor() = default;
		
		void open_preprocessing_result_file(char const *path);
		void process_and_output();
		void variant_graph_generator_will_handle_subgraph(libbio::variant const &, std::size_t const, std::size_t const);
	};
}


namespace {
	
	std::size_t progress_indicator_delegate::progress_step_max() const
	{
		return m_processor->m_preprocessing_result.handled_line_numbers.size();
	}
	
	
	std::size_t progress_indicator_delegate::progress_current_step() const
	{
		return m_processor->m_generator.processed_count();
	}
	
	
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
}


namespace vcf2multialign {
	
	void variant_graph_processor::open_preprocessing_result_file(char const *path)
	{
		lb::file_istream input_preprocessing_result_stream;
		lb::open_file_for_reading(path, input_preprocessing_result_stream);
		cereal::PortableBinaryInputArchive iarchive(input_preprocessing_result_stream);
		iarchive(m_preprocessing_result);
	}
	
	
	void variant_graph_processor::variant_graph_generator_will_handle_subgraph(
		lb::variant const &first_var,
		std::size_t const variant_count,
		std::size_t const path_count
	)
	{
#if 0
		auto const var_pos(first_var.zero_based_pos());
		dispatch_async_main(^{
			std::cerr << "\33[1G" << "Subgraph start: " << var_pos << " variants: " << variant_count << " paths: " << path_count << '\n';
		});
#endif
	}
	
	
	void variant_graph_processor::process_and_output()
	{
		try
		{
			m_generator = variant_graph_generator(*this, this->m_vcf_reader, this->m_reference, m_preprocessing_result);
			
			// Partition.
			progress_indicator_delegate indicator_delegate(*this);
			this->install_progress_indicator();
			
			dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Creating the variant graph…\n"; });
			this->progress_indicator().log_with_progress_bar("\t", indicator_delegate);
			m_generator.generate_graph();
			this->end_logging();
			this->uninstall_progress_indicator();
			
			dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Done.\n"; }); // FIXME: log statistics?
			
			// Output.
			cereal::PortableBinaryOutputArchive archive(this->m_output_stream);
			archive(m_generator.variant_graph());
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
	
	
	void create_variant_graph(
		char const *reference_file_path,
		char const *variant_file_path,
		char const *preprocessing_result_file_path,
		char const *output_graph_path,
		char const *reference_seq_name,
		bool const should_overwrite_files
	)
	{
		// main() calls dispatch_main() if this function does not throw or call std::exit or something similar.
		try
		{
			// Since the processor has Boost’s streams, it cannot be moved. Hence the use of a pointer.
			auto processor_ptr(std::make_unique <variant_graph_processor>());
			auto &processor(*processor_ptr);
			
			// These will eventually call std::exit if the file in question cannot be opened.
			processor.open_variants_file(variant_file_path);
			processor.open_preprocessing_result_file(preprocessing_result_file_path);
			processor.read_reference(reference_file_path, reference_seq_name);
			processor.open_output_file(output_graph_path, should_overwrite_files);
			
			processor.prepare_reader();
			
			// Run in background in order to be able to update a progress bar.
			lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
				[
					processor_ptr = std::move(processor_ptr)
				](){
					processor_ptr->process_and_output();
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
