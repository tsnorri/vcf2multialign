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
	
	
	class preprocess_logger : public v2m::variant_partitioner_delegate
	{
	protected:
		lb::dispatch_ptr <dispatch_queue_t>				m_queue;
		lb::file_ostream								m_log_output_stream;
		std::set <std::pair <std::size_t, std::size_t>>	m_reported_overlaps;
		std::mutex										m_mutex{};
		
	public:
		void open_log_file(char const *path, bool const should_overwrite_files);
		
		void variant_processor_no_field_for_identifier(std::string const &identifier) override;
		void variant_processor_found_variant_with_position_greater_than_reference_length(libbio::transient_variant const &var) override;
		void variant_processor_found_variant_with_no_suitable_alts(libbio::transient_variant const &var) override;
		void variant_processor_found_filtered_variant(libbio::transient_variant const &var, libbio::vcf_info_field_base const &field) override;
		void variant_processor_found_variant_with_ref_mismatch(libbio::transient_variant const &var, std::string_view const &ref_sub) override;
		void sample_sorter_found_overlapping_variant(lb::variant const &var, std::size_t const sample_idx, std::size_t const prev_end_pos) override;
		
	protected:
		template <typename t_extra>
		void log_wt(std::size_t const lineno, std::size_t const sample_idx, char const *reason, t_extra const &extra);
		
		template <typename t_extra>
		void log_wt(std::size_t const lineno, char const *reason, t_extra const &extra);
		
		void log_wt(std::size_t const lineno, std::size_t const sample_idx, char const *reason) { log_wt(lineno, sample_idx, reason, '.'); }
		void log_wt(std::size_t const lineno, char const *reason) { log_wt(lineno, reason, '.'); }
	};
	
	
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
	
	
	void preprocess_logger::open_log_file(char const *output_path, bool const should_overwrite_files)
	{
		auto const mode(lb::make_writing_open_mode({
			lb::writing_open_mode::CREATE,
			(should_overwrite_files ? lb::writing_open_mode::OVERWRITE : lb::writing_open_mode::NONE)
		}));
		lb::open_file_for_writing(output_path, m_log_output_stream, mode);
		
		m_queue.reset(dispatch_queue_create("fi.iki.tsnorri.vcf2multialign.preprocess_logger", DISPATCH_QUEUE_SERIAL));
		dispatch_async(*m_queue, ^{
			this->m_log_output_stream << "LINENO\tSAMPLE\tREASON\tEXTRA\n";
		});
	}
	
	
	void preprocess_logger::variant_processor_no_field_for_identifier(std::string const &identifier)
	{
		lb::dispatch_async_fn(dispatch_get_main_queue(), [identifier](){
			std::cerr << "\nWARNING: Did not find a field for identifier “" << identifier << "”.\n";
		});
	}
	
	
	void preprocess_logger::variant_processor_found_variant_with_position_greater_than_reference_length(libbio::transient_variant const &var)
	{
		auto const lineno(var.lineno());
		dispatch_async(dispatch_get_main_queue(), ^{
			std::cerr << "\nERROR: Found a variant with a position greater than the reference length on line " << lineno << "\n";
		});
	}
	
	
	void preprocess_logger::variant_processor_found_variant_with_no_suitable_alts(libbio::transient_variant const &var)
	{
		if (!m_queue)
			return;
		
		auto const lineno(var.lineno());
		dispatch_async(*m_queue, ^{
			this->log_wt(lineno, "Variant has no ALTs that could be handled");
		});
	}
	
	
	void preprocess_logger::variant_processor_found_filtered_variant(libbio::transient_variant const &var, libbio::vcf_info_field_base const &field)
	{
		if (!m_queue)
			return;
		
		auto const lineno(var.lineno());
		auto const &field_id(field.get_metadata()->get_id());
		dispatch_async(*m_queue, ^{
			this->log_wt(lineno, "Variant has a filtered field set", field_id);
		});
	}
	
	
	void preprocess_logger::variant_processor_found_variant_with_ref_mismatch(libbio::transient_variant const &var, std::string_view const &ref_sub)
	{
		auto const lineno(var.lineno());
		std::string var_ref(var.ref());
		std::string ref_sub_copy(ref_sub);
		lb::dispatch_async_fn(dispatch_get_main_queue(), [lineno, var_ref = std::move(var_ref), ref_sub_copy = std::move(ref_sub_copy)](){
			std::cerr << "\nWARNING: reference column mismatch on line " << lineno << ": expected '" << ref_sub_copy << "', got '" << var_ref<< "'\n";
		});
	}
	
	
	void preprocess_logger::sample_sorter_found_overlapping_variant(lb::variant const &var, std::size_t const sample_idx, std::size_t const prev_end_pos)
	{
		if (!m_queue)
			return;
		
		auto const lineno(var.lineno());
		auto const p(std::make_pair(lineno, sample_idx));
		bool should_write(false);
		
		{
			std::lock_guard lock(m_mutex);
			if (this->m_reported_overlaps.end() == this->m_reported_overlaps.find(p))
				should_write = true;
		}
		
		if (should_write)
		{
			dispatch_async(*m_queue, ^{
				this->m_reported_overlaps.insert(p);
				this->log_wt(lineno, sample_idx, "Overlapping alternatives", prev_end_pos);
			});
		}
	}
	
	
	template <typename t_extra>
	void preprocess_logger::log_wt(std::size_t const lineno, std::size_t const sample_idx, char const *reason, t_extra const &extra)
	{
		this->m_log_output_stream << lineno << '\t' << sample_idx << '\t' << reason << '\t' << extra << '\n';
	}
	
	
	template <typename t_extra>
	void preprocess_logger::log_wt(std::size_t const lineno, char const *reason, t_extra const &extra)
	{
		this->m_log_output_stream << lineno << "\t.\t" << reason << '\t' << extra << '\n';
	}
}


namespace vcf2multialign {
	
	class cut_position_processor final :	public vcf_processor,
											public output_stream_controller,
											public reference_controller,
											public progress_indicator_manager,
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
			variant_partitioner partitioner(
				*this,
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
			auto processor_ptr(std::make_unique <cut_position_processor>(std::move(field_names_for_filter_if_set), chr_name, minimum_subgraph_distance));
			auto &processor(*processor_ptr);
			
			// These will eventually call std::exit if the file in question cannot be opened.
			processor.open_variants_file(variant_file_path);
			processor.read_reference(reference_path, reference_seq_name);
			processor.open_output_file(output_path, should_overwrite_files);
			if (log_path)
				processor.open_log_file(log_path, should_overwrite_files);
			
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
