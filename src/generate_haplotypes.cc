/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/map.hpp>
#include <iostream>
#include <thread>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/generate_haplotypes.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/tasks/all_haplotypes_task.hh>
#include <vcf2multialign/tasks/preparation_task.hh>
#include <vcf2multialign/tasks/reduce_samples_task.hh>

namespace ios	= boost::iostreams;
namespace v2m	= vcf2multialign;


namespace {
	
	class generate_context final :
		public v2m::all_haplotypes_task_delegate,
		public v2m::preparation_task_delegate,
		public v2m::reduce_samples_task_delegate
	{
	protected:
		typedef std::set <
			std::unique_ptr <v2m::task>,
			v2m::pointer_cmp <v2m::task>
		> task_set;
		
	protected:
		std::mutex											m_tasks_mutex{};
		task_set											m_tasks;
		v2m::vector_type									m_reference;
		v2m::vcf_mmap_input									m_vcf_input;
		v2m::mmap_handle									m_vcf_handle;
		
		v2m::error_logger									m_error_logger;
		v2m::status_logger									m_status_logger;
		
		v2m::ploidy_map										m_ploidy;
		v2m::variant_set									m_skipped_variants;
		v2m::subgraph_map									m_subgraph_starting_points;	// Variant line numbers by character pointers.
		v2m::alt_checker									m_alt_checker;
		
		boost::optional <std::string>						m_out_reference_fname;
		std::string											m_null_allele_seq;
		v2m::sv_handling									m_sv_handling_method;
		std::size_t											m_chunk_size{0};
		std::size_t											m_min_path_length{0};
		std::size_t											m_generated_path_count{0};
		
		bool												m_should_overwrite_files{false};
		bool												m_should_reduce_samples{false};
		
	public:
		generate_context(
			char const *out_reference_fname,
			char const *null_allele_seq,
			v2m::sv_handling const sv_handling_method,
			std::size_t const chunk_size,
			std::size_t const min_path_length,
			std::size_t const generated_path_count,
			bool const should_overwrite_files,
			bool const should_reduce_samples
		):
			m_vcf_input(m_vcf_handle),
			m_null_allele_seq(null_allele_seq),
			m_sv_handling_method(sv_handling_method),
			m_chunk_size(chunk_size),
			m_min_path_length(min_path_length),
			m_generated_path_count(generated_path_count),
			m_should_overwrite_files(should_overwrite_files),
			m_should_reduce_samples(should_reduce_samples)
		{
			finish_init(out_reference_fname);
		}
	
		generate_context(generate_context const &) = delete;
		generate_context(generate_context &&) = delete;
		
		void cleanup() { m_status_logger.uninstall(); delete this; }
		
		void load_and_generate(
			char const *reference_fname,
			char const *variants_fname,
			char const *report_fname,
			bool const should_check_ref
		);
			
		void store_and_execute_in_group(std::unique_ptr <v2m::task> &&task, dispatch_group_t group);
		
		// preparation_task_delegate
		virtual void task_did_finish(v2m::preparation_task &task) override;
		
		// all_haplotypes_task_delegate
		virtual void task_did_finish(v2m::all_haplotypes_task &task) override;
		
		// reduce_samples_task_delegate
		virtual void task_did_finish(v2m::reduce_samples_task &task) override;
		virtual void store_and_execute(std::unique_ptr <v2m::task> &&task) override;
		
	protected:
		void finish_init(char const *out_reference_fname);
		v2m::task *store_task(std::unique_ptr <v2m::task> &&task);
		void remove(v2m::task &task);
	};
	
	
	void generate_context::finish_init(
		char const *out_reference_fname
	)
	{
		if (out_reference_fname)
			m_out_reference_fname.emplace(out_reference_fname);
	}
	
	
	v2m::task *generate_context::store_task(std::unique_ptr <v2m::task> &&task)
	{
		std::lock_guard <std::mutex> lock_guard(m_tasks_mutex);
		auto st(m_tasks.emplace(std::move(task)));
		
		// task is now invalid.
		v2m::always_assert(st.second);
		v2m::task *inserted_task(st.first->get());
		
		return inserted_task;
	}
	
	
	void generate_context::store_and_execute(std::unique_ptr <v2m::task> &&task)
	{
		auto inserted_task(store_task(std::move(task)));
		auto queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
		v2m::dispatch(inserted_task).async <&v2m::task::execute>(queue);
	}
	
	
	void generate_context::store_and_execute_in_group(std::unique_ptr <v2m::task> &&task, dispatch_group_t group)
	{
		auto inserted_task(store_task(std::move(task)));
		auto queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
		v2m::dispatch(inserted_task).group_async <&v2m::task::execute>(group, queue);
	}
	
	
	void generate_context::remove(v2m::task &task)
	{
		// Use the C++14 templated find.
		auto it(m_tasks.find(task));
		v2m::always_assert(m_tasks.cend() != it);
		m_tasks.erase(it);
	}
	
	
	void generate_context::load_and_generate(
		char const *reference_fname,
		char const *variants_fname,
		char const *report_fname,
		bool const should_check_ref
	)
	{
		m_status_logger.install();
		
		// It is easier to open the files here rather than in the preparation task.
		v2m::vcf_reader reader(m_vcf_input);
		
		{
			v2m::file_istream ref_fasta_stream;
			
			v2m::open_file_for_reading(reference_fname, ref_fasta_stream);
			m_vcf_handle.open(variants_fname);
			m_vcf_input.reset_range();
			
			if (report_fname)
			{
				m_error_logger.prepare();
				v2m::open_file_for_writing(
					report_fname,
					m_error_logger.output_stream(),
					m_should_overwrite_files
				);
				m_error_logger.write_header();
			}
			
			std::cerr << "Reading the VCF header…" << std::endl;
			reader.read_header();
			
			// Read the reference file and place its contents into reference.
			// If minimu path length was not given, set it to the square root of the reference sequence length.
			v2m::read_single_fasta_seq(ref_fasta_stream, m_reference);
			if (m_should_reduce_samples && 0 == m_min_path_length)
			{
				m_min_path_length = std::ceil(std::sqrt(m_reference.size()));
				std::cerr << "Set minimum path length to " << m_min_path_length << '.' << std::endl;
			}
		}
		
		std::unique_ptr <v2m::task> task(
			new v2m::preparation_task(
				*this,
				m_status_logger,
				m_error_logger,
				m_reference,
				std::move(reader),
				m_sv_handling_method,
				should_check_ref
			)
		);
		
		store_and_execute(std::move(task));
	}
	
	
	void generate_context::task_did_finish(v2m::preparation_task &task)
	{
		v2m::vcf_reader reader(std::move(task.vcf_reader()));
		
		m_ploidy = std::move(task.ploidy_map());
		m_skipped_variants = std::move(task.skipped_variants());
		m_alt_checker = std::move(task.alt_checker());
		m_subgraph_starting_points = std::move(task.subgraph_starting_points());
		auto const record_count(task.step_count());
		auto const records_with_valid_alts(m_alt_checker.records_with_valid_alts());
		
		remove(task);
		// task is now invalid.
		
		if (m_should_reduce_samples)
		{
			auto const sample_ploidy_sum(boost::accumulate(m_ploidy | boost::adaptors::map_values, 0));
			auto const hw_concurrency(std::thread::hardware_concurrency());
			
			std::unique_ptr <v2m::task> task(
				new v2m::reduce_samples_task(
					*this,
					m_status_logger,
					m_error_logger,
					hw_concurrency,
					std::move(reader),
					m_vcf_handle,
					m_reference,
					m_null_allele_seq,
					m_alt_checker,
					m_subgraph_starting_points,
					m_skipped_variants,
					m_out_reference_fname,
					m_sv_handling_method,
					record_count,
					m_generated_path_count,
					sample_ploidy_sum,
					m_min_path_length,
					m_should_overwrite_files
				)
			);
			
			store_and_execute(std::move(task));
		}
		else
		{
			// FIXME: create the queues in a function? The case above (m_should_reduce_samples) is
			// going to be handled with a custom variant handler (with enumerate_sample_genotypes
			// replaced but everything else should be OK as is) and with a custom task.
			auto concurrent_queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
			v2m::dispatch_ptr <dispatch_queue_t> parsing_queue(concurrent_queue, true);
			v2m::dispatch_ptr <dispatch_queue_t> worker_queue(
				dispatch_queue_create("fi.iki.tsnorri.vcf2multialign.worker_queue", DISPATCH_QUEUE_SERIAL),
				false
			);
			
			std::unique_ptr <v2m::task> task(
				new v2m::all_haplotypes_task(
					*this,
					worker_queue,
					m_status_logger,
					m_error_logger,
					std::move(reader),
					m_alt_checker,
					m_reference,
					m_null_allele_seq,
					m_ploidy,
					m_skipped_variants,
					m_out_reference_fname,
					m_sv_handling_method,
					record_count,
					records_with_valid_alts,
					m_chunk_size,
					m_should_overwrite_files
				)
			);
			
			store_and_execute(std::move(task));
		}
	}
	
	
	void generate_context::task_did_finish(v2m::reduce_samples_task &task)
	{
		cleanup();
		exit(EXIT_SUCCESS);
	}
	
	
	void generate_context::task_did_finish(v2m::all_haplotypes_task &task)
	{
		// After calling cleanup *this is no longer valid.
		//std::cerr << "Calling cleanup" << std::endl;
		cleanup();
	
		//std::cerr << "Calling exit" << std::endl;
		exit(EXIT_SUCCESS);
	}
}


namespace vcf2multialign {
	
	void generate_haplotypes(
		char const *reference_fname,
		char const *variants_fname,
		char const *out_reference_fname,
		char const *report_fname,
		char const *null_allele_seq,
		std::size_t const chunk_size,
		std::size_t const min_path_length,
		std::size_t const generated_path_count,
		sv_handling const sv_handling_method,
		bool const should_overwrite_files,
		bool const should_check_ref,
		bool const should_reduce_samples
	)
	{
		// generate_context needs to be allocated on the heap because later dispatch_main is called.
		// The class deallocates itself in cleanup().
		generate_context *ctx(new generate_context(
			out_reference_fname,
			null_allele_seq,
			sv_handling_method,
			chunk_size,
			min_path_length,
			generated_path_count,
			should_overwrite_files,
			should_reduce_samples
		));
			
		ctx->load_and_generate(
			reference_fname,
			variants_fname,
			report_fname,
			should_check_ref
		);
	}
}
