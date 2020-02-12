/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cereal/archives/portable_binary.hpp>
#include <libbio/progress_indicator.hh>
#include <vcf2multialign/graph/variant_graph_generator.hh>
#include <vcf2multialign/preprocess/preprocess_logger.hh>
#include <vcf2multialign/preprocess/variant_partitioner.hh>
#include <vcf2multialign/utility/check_ploidy.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/utility/progress_indicator_manager.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "create_variant_graph.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	class variant_graph_context; // Fwd.
}


namespace {
	
	class progress_indicator_delegate final : public lb::progress_indicator_delegate
	{
	protected:
		v2m::variant_graph_generator const * const	m_generator{};
		std::size_t									m_progress_step_max{};
		
	public:
		progress_indicator_delegate(
			v2m::variant_graph_generator const &generator,
			std::size_t const progress_step_max
		):
			m_generator(&generator),
			m_progress_step_max(progress_step_max)
		{
		}
		
		virtual void progress_log_extra() const override {};
		inline virtual std::size_t progress_step_max() const override { return m_progress_step_max; }
		inline virtual std::size_t progress_current_step() const override { return m_generator->processed_count(); }
	};
}


namespace vcf2multialign {
	
	class variant_graph_context :	public vcf_processor,
									public output_stream_controller,
									public reference_controller,
									public progress_indicator_manager,
									public virtual variant_graph_generator_delegate
	{
		friend progress_indicator_delegate;
		
	protected:
		preprocessing_result		m_preprocessing_result;	// Needed for variant_graph_generator.
		
	public:
		variant_graph_context() = default;
		
		virtual class variant_graph_generator &variant_graph_generator() = 0;
		void process_and_output();
		void variant_graph_generator_will_handle_subgraph(libbio::variant const &, std::size_t const, std::size_t const);
		
	protected:
		virtual std::size_t progress_step_max() const = 0;
		virtual void generate_variant_graph(progress_indicator_delegate &indicator_delegate) = 0;
		virtual void log_statistics() const = 0;
	};
	
	
	class variant_graph_precalculated_context final : public variant_graph_context
	{
	protected:
		variant_graph_precalculated_generator	m_generator;
		
	public:
		variant_graph_precalculated_context() = default;
		
		class variant_graph_precalculated_generator &variant_graph_generator() override { return m_generator; }
		void open_preprocessing_result_file(char const *path);
		void prepare_generator();
		
	protected:
		std::size_t progress_step_max() const override { return m_preprocessing_result.handled_line_numbers.size(); }
		void generate_variant_graph(progress_indicator_delegate &indicator_delegate) override;
		virtual void log_statistics() const override;
	};
	
	
	class variant_graph_single_pass_context final : public variant_graph_context,
													public virtual variant_graph_single_pass_generator_delegate,
													public preprocess_logger // For the case in which preprocessing was not done.
													{
	protected:
		std::string								m_chromosome_name;
		std::vector <std::string>				m_field_names_for_filter_if_set;
		variant_graph_single_pass_generator		m_generator;
		
	public:
		variant_graph_single_pass_context(std::string &&chromosome_name, std::vector <std::string> &&field_names_for_filter_if_set):
			m_chromosome_name(std::move(chromosome_name)),
			m_field_names_for_filter_if_set(std::move(field_names_for_filter_if_set))
		{
		}
		
		class variant_graph_single_pass_generator &variant_graph_generator() override { return m_generator; }
		std::tuple <std::size_t, std::uint8_t> check_ploidy_from_vcf();
		void prepare_generator();
		
	protected:
		std::size_t progress_step_max() const override { return 0; }
		void generate_variant_graph(progress_indicator_delegate &indicator_delegate) override;
		virtual void log_statistics() const override;
	};
}


namespace {
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
}


namespace vcf2multialign {
	
	void variant_graph_context::variant_graph_generator_will_handle_subgraph(
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
	
	
	void variant_graph_precalculated_context::open_preprocessing_result_file(char const *path)
	{
		lb::file_istream input_preprocessing_result_stream;
		lb::open_file_for_reading(path, input_preprocessing_result_stream);
		cereal::PortableBinaryInputArchive iarchive(input_preprocessing_result_stream);
		iarchive(m_preprocessing_result);
	}
	
	
	std::tuple <std::size_t, std::uint8_t> variant_graph_single_pass_context::check_ploidy_from_vcf()
	{
		// Check the ploidy.
		ploidy_map ploidy;
		check_ploidy(this->vcf_reader(), ploidy);
		auto const donor_count(ploidy.size());
		if (!donor_count)
		{
			std::cerr << "WARNING: No donors found." << std::endl;
			std::exit(EXIT_SUCCESS);
		}
		auto const chr_count(ploidy.begin()->second);
		
		return std::make_tuple(donor_count, chr_count);
	}
	
	
	void variant_graph_precalculated_context::prepare_generator()
	{
		m_generator = variant_graph_precalculated_generator(*this, this->m_vcf_reader, this->m_reference, m_preprocessing_result);
		m_generator.sample_sorter().set_sample_indexer(m_generator.sample_indexer()); // Fix the pointer. FIXME: should be done in variant_graph_generator’s move assignment operator.
	}
	
	
	void variant_graph_single_pass_context::prepare_generator()
	{
		auto const [donor_count, chr_count] = check_ploidy_from_vcf();
		m_generator = variant_graph_single_pass_generator(
			*this,
			this->m_vcf_reader,
			this->m_reference,
			m_chromosome_name,
			donor_count,
			chr_count
		);
		m_generator.sample_sorter().set_sample_indexer(m_generator.sample_indexer()); // Fix the pointer. FIXME: should be done in variant_graph_generator’s move assignment operator.
	}
	
	
	void variant_graph_precalculated_context::generate_variant_graph(progress_indicator_delegate &indicator_delegate)
	{
		dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Creating the variant graph…\n"; });
		this->progress_indicator().log_with_progress_bar("\t", indicator_delegate);
		m_generator.generate_graph();
	}
	
	
	void variant_graph_single_pass_context::generate_variant_graph(progress_indicator_delegate &indicator_delegate)
	{
		this->progress_indicator().log_with_counter(lb::copy_time() + "Creating the variant graph…", indicator_delegate);
		m_generator.generate_graph(m_field_names_for_filter_if_set);
	}


	void variant_graph_precalculated_context::log_statistics() const
	{
		dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Done.\n"; }); // FIXME: log statistics?
	}


	void variant_graph_single_pass_context::log_statistics() const
	{
		dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Done. Chromosome ID mismatches: " << this->chrom_id_mismatches() << "\n"; }); // FIXME: log more statistics?
	}
	
	
	void variant_graph_context::process_and_output()
	{
		try
		{
			// Partition.
			progress_indicator_delegate indicator_delegate(this->variant_graph_generator(), this->progress_step_max());
			this->install_progress_indicator();
			this->generate_variant_graph(indicator_delegate);
			this->end_logging();
			this->uninstall_progress_indicator();
			this->log_statistics();
			
			// Output.
			cereal::PortableBinaryOutputArchive archive(this->m_output_stream);
			archive(this->variant_graph_generator().variant_graph());
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
	
	
	void create_variant_graph_preprocessed(
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
			// Since the context has Boost’s streams, it cannot be moved. Hence the use of a pointer.
			auto ctx_ptr(std::make_unique <variant_graph_precalculated_context>());
			auto &ctx(*ctx_ptr);
			
			ctx.open_variants_file(variant_file_path);
			ctx.open_preprocessing_result_file(preprocessing_result_file_path);
			ctx.read_reference(reference_file_path, reference_seq_name);
			ctx.open_output_file(output_graph_path, should_overwrite_files);
			ctx.prepare_reader();
			ctx.prepare_generator();
		
			// Run in background in order to be able to update a progress bar.
			lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
				[ctx_ptr = std::move(ctx_ptr)](){
					ctx_ptr->process_and_output();
				}
			);
		}
		catch (lb::assertion_failure_exception const &exc)
		{
			log_assertion_failure_exception(exc);
			throw exc;
		}
	}
	
	
	void create_variant_graph_single_pass(
		char const *reference_file_path,
		char const *variant_file_path,
		char const *output_graph_path,
		char const *log_path,
		char const *reference_seq_name,
		char const *chr_id,
		std::vector <std::string> &&field_names_for_filter_if_set,
		bool const should_overwrite_files
	)
	{
		// main() calls dispatch_main() if this function does not throw or call std::exit or something similar.
		try
		{
			// Since the context has Boost’s streams, it cannot be moved. Hence the use of a pointer.
			auto ctx_ptr(std::make_unique <variant_graph_single_pass_context>(chr_id, std::move(field_names_for_filter_if_set)));
			auto &ctx(*ctx_ptr);
			
			ctx.open_variants_file(variant_file_path);
			ctx.read_reference(reference_file_path, reference_seq_name);
			ctx.open_output_file(output_graph_path, should_overwrite_files);
			ctx.open_log_file(log_path, should_overwrite_files);
			ctx.prepare_reader();
			ctx.prepare_generator();
			
			// Run in background in order to be able to update a progress bar.
			lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
				[ctx_ptr = std::move(ctx_ptr)](){
					ctx_ptr->process_and_output();
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
