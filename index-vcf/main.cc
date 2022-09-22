/*
 * Copyright (c) 2019–2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdlib>
#include <cereal/archives/portable_binary.hpp>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/dispatch.hh>
#include <libbio/progress_indicator.hh>
#include <unistd.h>
#include <vcf2multialign/utility/dispatch_cli_runner.hh>
#include <vcf2multialign/utility/dispatch_exit_guard.hh>
#include <vcf2multialign/utility/log_assertion_failure.hh>
#include <vcf2multialign/vcf_index.hh>
#include <vcf2multialign/vcf_processor.hh>
#include "cmdline.h"


namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace {
	
	class vcf_indexer; // Fwd.
	
	
	class progress_indicator_delegate final : public lb::progress_indicator_delegate
	{
	protected:
		vcf_indexer const * const	m_indexer{};
		
	public:
		progress_indicator_delegate(vcf_indexer const &indexer):
			m_indexer(&indexer)
		{
		}
		
		virtual std::size_t progress_step_max() const override { return 0; }
		virtual std::size_t progress_current_step() const override;
		virtual void progress_log_extra() const override {};
	};
	
	
	class vcf_indexer :	public v2m::dispatch_cli_runner,
						public v2m::vcf_processor,
						public v2m::output_stream_controller
	{
	protected:
		std::string		m_chromosome_name;
		v2m::vcf_index	m_index;
		
	public:
		vcf_indexer() = default;
		
		vcf_indexer(char const *m_chromosome_name):
			m_chromosome_name(m_chromosome_name ?: "")
		{
		}
		
		std::size_t processed_count() const { return this->m_vcf_reader.counter_value(); }
		
	protected:
		void do_work() override;
		void generate_index(progress_indicator_delegate &progress_delegate);
	};
	
	
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
	
	
	std::size_t progress_indicator_delegate::progress_current_step() const
	{
		return m_indexer->processed_count();
	}
	
	
	void vcf_indexer::generate_index(progress_indicator_delegate &progress_delegate)
	{
		bool should_continue(false);
		auto &reader(this->m_vcf_reader);
		reader.set_parsed_fields(vcf::field::POS);
		v2m::vcf_index::chromosome_entry *entry_ptr{};
		std::size_t prev_pos{};
		this->m_vcf_reader.parse(
			[
				this,
				&entry_ptr,
				&prev_pos
			](vcf::transient_variant const &var) -> bool
			{
				auto const &chrom_id(var.chrom_id());
				if (!m_chromosome_name.empty() && m_chromosome_name != chrom_id)
					return true;
				
				auto const var_pos(var.pos());
				auto const var_offset(m_vcf_reader.variant_offset());
				
				// Check if the chromosome identifier has changed. If not,
				// check if the position changed in order to store only the first entry for
				// each position.
				if (!entry_ptr || chrom_id != entry_ptr->chromosome_identifier)
				{
					entry_ptr = &m_index.add_entry(chrom_id);
					entry_ptr->add_position(var_pos - 1, var_offset);
				}
				else if (prev_pos < var_pos)
				{
					entry_ptr->add_position(var_pos - 1, var_offset);
				}
				
				prev_pos = var_pos;
				return true;
			}
		);
	}
	
	
	void vcf_indexer::do_work()
	{
		progress_indicator_delegate indicator_delegate(*this);
		this->progress_indicator().log_with_counter(lb::copy_time() + "Indexing the variant file…", indicator_delegate);
		this->generate_index(indicator_delegate);
		this->end_logging();
		this->uninstall_progress_indicator();
		
		dispatch_async_main(^{ lb::log_time(std::cerr); std::cerr << "Done.\n"; }); // FIXME: log statistics?
		
		// Output.
		cereal::PortableBinaryOutputArchive archive(this->m_output_stream);
		archive(m_index);
		this->m_output_stream << std::flush;
	}
	
	
	void index_vcf(
		char const *variant_file_path,
		char const *output_index_path,
		char const *chromosome_name,
		bool const should_overwrite_files
	)
	{
		typedef v2m::dispatch_exit_guard_helper <vcf_indexer> wrapped_indexer_type;
		auto indexer_ptr(std::make_unique <wrapped_indexer_type>(chromosome_name));
		auto &indexer(indexer_ptr->value);
		indexer.open_variants_file(variant_file_path);
		indexer.open_output_file(output_index_path, should_overwrite_files);
		indexer.prepare_reader();
		
		// Run in background in order to be able to update a progress bar.
		lb::dispatch_async_fn(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0),
			[indexer_ptr = std::move(indexer_ptr)](){
				indexer_ptr->value.run(true);
			}
		);
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
	
	try
	{
		index_vcf(args_info.variants_arg, args_info.output_arg, args_info.chromosome_arg, args_info.overwrite_flag);
	}
	catch (lb::assertion_failure_exception const &exc)
	{
		v2m::log_assertion_failure_exception(exc);
	}
	
	dispatch_main();
	// Not reached.
	return EXIT_SUCCESS;
}
