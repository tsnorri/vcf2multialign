/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/io/ios_state.hpp>
#include <vcf2multialign/file_handling.hh>
#include <vcf2multialign/tasks/read_subgraph_variants_task.hh>
#include <vcf2multialign/tasks/reduce_samples_task.hh>


namespace vcf2multialign { namespace detail {
	void reduce_samples_progress_counter::calculate_step_count(
		std::size_t const valid_record_count,
		std::size_t const path_count,
		std::size_t const merge_tasks
	)
	{
		std::size_t step_count(0);
	
		// Handling each variant in read_subgraph_variants_task.
		// FIXME: this does not take into account generating the paths.
		m_rsv_steps = valid_record_count;
		m_step_count += m_rsv_steps;
	
		// Calculating the weight of each edge between subgraphs.
		m_edge_weight_steps = merge_tasks * path_count * path_count;
		m_step_count += m_edge_weight_steps;
		
		// Merging.
		m_merge_tasks = merge_tasks;
		m_step_count += m_merge_tasks;
	}
}}


namespace vcf2multialign {
	
	void reduce_samples_task::init_read_subgraph_variants_task(
		read_subgraph_variants_task &task,
		std::size_t const task_idx,
		graph_range &&range
	)
	{
		auto const qname(boost::str(boost::format("fi.iki.tsnorri.vcf2multialign.worker-queue-%u") % task_idx));
		dispatch_ptr <dispatch_queue_t> worker_queue(
			dispatch_queue_create(qname.c_str(), DISPATCH_QUEUE_SERIAL),
			false
		);
			
		read_subgraph_variants_task temp_task(
			*this,
			worker_queue,
			*m_status_logger,
			*m_error_logger,
			m_reader,
			*m_vcf_input_handle,
			*m_alt_checker,
			*m_reference,
			m_sv_handling_method,
			*m_skipped_variants,
			*m_out_reference_fname,
			std::move(range),
			task_idx,
			m_generated_path_count,
			m_sample_ploidy_sum
		);
		
		task = std::move(temp_task);
	}
	
	
	void reduce_samples_task::start_merge_task(std::size_t const lhs_idx, std::size_t const rhs_idx)
	{
		auto const &lhs(m_subgraphs[lhs_idx]);
		auto const &rhs(m_subgraphs[rhs_idx]);
		std::unique_ptr <task> task(new merge_subgraph_paths_task(*this, *m_status_logger, lhs, rhs, lhs_idx, m_generated_path_count));
		
		m_delegate->store_and_execute(std::move(task));
	}
	
	
	auto reduce_samples_task::create_haplotype(
		std::size_t const sample_no,
		std::size_t const ploidy
	) -> haplotype_map_type::iterator
	{
		return m_haplotypes.emplace(
			std::piecewise_construct,
			std::forward_as_tuple(sample_no),
			std::forward_as_tuple(ploidy)
		).first;
	}
	
	
	void reduce_samples_task::prepare_haplotypes()
	{
		// FIXME: partially duplicate code with all_haplotypes_task.
		// FIXME: m_generated_path_count (plus some small value) may not exceed the maximum number of open files.
		for (std::size_t i(0); i < m_generated_path_count; ++i)
		{
			auto const fname(boost::str(boost::format("%u") % (1 + i)));
			auto it(create_haplotype(1 + i, 1));
			auto &haplotype_vec(it->second);
			open_file_channel_for_writing(
				fname.c_str(),
				haplotype_vec[0].output_stream,
				m_write_semaphore,
				m_should_overwrite_files
			);
		}
		
		if (m_out_reference_fname->operator bool())
		{
			always_assert(
				m_haplotypes.cend() == m_haplotypes.find(REF_SAMPLE_NUMBER),
				"REF_SAMPLE_NUMBER already in use"
			);
			
			auto it(create_haplotype(REF_SAMPLE_NUMBER, 1));
			auto &haplotype_vec(it->second);
			open_file_channel_for_writing(
				(*m_out_reference_fname)->c_str(),
				haplotype_vec[0].output_stream,
				m_write_semaphore,
				m_should_overwrite_files
			);
		}
	}
	
	
	void reduce_samples_task::task_did_finish(
		merge_subgraph_paths_task &task,
		std::vector <reduced_subgraph::path_index> &&matchings,
		merge_subgraph_paths_task::weight_type const matching_weight
	)
	{
		m_progress_counter.merge_subgraph_paths_task_did_finish();
		m_matching_weight += matching_weight;
		
		// As per [container.requirements.dataraces] this should be thread-safe since different
		// threads may not modify the same element.
		auto const idx(task.left_subgraph_index());
		m_path_matchings[idx] = std::move(matchings);
		
		// According to CppReference, the above should become a visible side effect of m_remaining_merge_tasks.load()
		// (a consequence of calling m_remaining_merge_tasks.operator--() and the default memory order being 
		// std::memory_order_seq_cst):
		// “–– everything that happened-before a store in one thread becomes a visible
		// side effect in the thread that did a load –– ”
		
		std::size_t const remaining_merge_tasks(--m_remaining_merge_tasks);
		m_status_logger->set_message(boost::str(boost::format("Reducing samples… (%d subgraphs to be merged)") % remaining_merge_tasks));
		
		if (0 == remaining_merge_tasks)
		{
			m_status_logger->finish_logging();
			m_progress_counter.reset_step_count(m_alt_checker->records_with_valid_alts());
			
			m_status_logger->log([weight = std::size_t(m_matching_weight)](){
				std::cerr << "Total matching weight was " << weight << '.' << std::endl;
			});
			
			m_status_logger->log_message_progress_bar("Writing the sequences…");
			auto task(new sequence_writer_task(
				*this,
				*m_status_logger,
				*m_error_logger,
				m_reader,
				*m_alt_checker,
				*m_reference,
				*m_null_allele_seq
			));
				
			std::unique_ptr <class task> task_ptr(task);
			
			prepare_haplotypes();
			task->prepare(m_haplotypes);
			
			m_subgraph_iterator = m_subgraphs.cbegin();
			
			m_delegate->store_and_execute(std::move(task_ptr));
		}
	}
	
	
	void reduce_samples_task::handled_all_haplotypes(sequence_writer_task &task)
	{
		m_status_logger->finish_logging();
	}
	
	
	void reduce_samples_task::task_did_finish(sequence_writer_task &task)
	{
		// Wait with a dispatch group that every file has finished writing.
		dispatch_queue_t queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
		dispatch_ptr <dispatch_group_t> group(dispatch_group_create());
		for (auto &kv : m_haplotypes)
		{
			// Actually there's only one element in each kv.second.
			for (auto &haplotype : kv.second)
			{
				// Make sure that dispatch_group_notify_fn does not fire before this
				// block has been executed.
				// FIXME: copying group may not be strictly needed since we call dispatch_group_async_fn with the same group.
				dispatch_group_async_fn(*group, queue, [&haplotype, group](){
					// Flush buffers and close. This is safe b.c.
					// set_closing_group() is called in the same thread as close().
					haplotype.output_stream->set_closing_group(group);
					haplotype.output_stream.close();
				});
			}
		}
		
		dispatch_group_notify_fn(*group, queue, [this](){
			// finish_logging() gets called in handled_all_haplotypes(). 
			m_delegate->task_did_finish(*this);
		});
	}
	
	
	void reduce_samples_task::enumerate_sample_genotypes(
		transient_variant const &var,
		std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb	// sample_no, chr_idx, alt_idx, is_phased
	)
	{
		auto const var_lineno(var.lineno());
		auto const first_var_lineno(1 + m_reader.last_header_lineno());
		
		// Find the correct subgraph.
		// Update m_path_permutation if the next subgraph is chosen.
		std::size_t subgraph_start_lineno(0);
		while (true)
		{
			always_assert(m_subgraphs.cend() != m_subgraph_iterator);
			
			auto const &subgraph(*m_subgraph_iterator);
			auto const subgraph_var_count(subgraph.variant_count());
			subgraph_start_lineno = subgraph.start_lineno();
			
			always_assert(subgraph_start_lineno <= var_lineno);
			if (subgraph.contains_var_lineno(var_lineno))
				break;
			
			// Update m_path_permutation.
			auto const subgraph_idx(m_subgraph_iterator - m_subgraphs.cbegin());
			auto const &matching(m_path_matchings[subgraph_idx]);
			auto const previous_permutation(m_path_permutation); // Copy.
			
			std::size_t i(0);
			for (auto const idx : previous_permutation)
			{
				m_path_permutation[i] = matching[idx];
				++i;
			}
			
			++m_subgraph_iterator;
		}
		
		// Enumerate the genotypes in the permutation order.
		auto const pos_in_subgraph(m_subgraph_iterator->seq_position(var_lineno));
		std::size_t sample_no(0);
		static_assert(0 == REF_SAMPLE_NUMBER);
		cb(REF_SAMPLE_NUMBER, 0, 0, true);	// REF.
		for (auto const path_idx : m_path_permutation)
		{
			auto const &sequence(m_subgraph_iterator->path_sequence(path_idx));
			auto const alt_idx(sequence[pos_in_subgraph]);
			
			// REF is zero, see the assertion above.
			cb(++sample_no, 0, alt_idx, true);
		}
		
		m_progress_counter.sequence_writer_task_did_handle_variant();
	}
	
	
	void reduce_samples_task::task_did_finish(read_subgraph_variants_task &task, reduced_subgraph &&rsg)
	{
		dispatch_semaphore_signal(*m_subgraph_semaphore);
		
		// Store the reduced subgraph. Use the mutex to make this thread-safe.
		std::lock_guard <std::mutex> lock_guard(m_subgraph_mutex);
		// Tasks are numbered in subgraph starting point order.
		auto const subgraph_idx(task.task_idx());
		m_subgraphs[subgraph_idx] = std::move(rsg);
		m_subgraph_bitmap[subgraph_idx] = true;
		
		// Start merging tasks if possible.
		// Currently start_merge_task requires the lock above.
		if (subgraph_idx && m_subgraph_bitmap[subgraph_idx - 1])
			start_merge_task(subgraph_idx - 1, subgraph_idx);
		
		auto const last_idx(m_subgraph_bitmap.size() - 1);
		if (subgraph_idx < last_idx && m_subgraph_bitmap[last_idx])
			start_merge_task(subgraph_idx, subgraph_idx + 1);
	}
	
	
	void reduce_samples_task::execute()
	{
		std::vector <graph_range> graph_ranges;
		
		// FIXME: the following block contains some duplicated code.
		// Create a task for each sufficiently long range.
		{
			std::vector <std::size_t> skipped_lines;
			auto const &input(m_reader.vcf_input());
			auto const buffer_start(input.buffer_start());
			auto const buffer_end(input.buffer_end());
			auto current_start(input.first_variant_start());
			std::size_t current_start_lineno(1 + m_reader.last_header_lineno());
			std::size_t current_lineno(current_start_lineno);	// Points to the last line that was checked for valid alts.
			std::size_t current_range_variant_count(0);
		
			for (auto const &kv : *m_subgraph_starting_points)
			{
				auto const next_subgraph_start(kv.first);
				auto const var_lineno(kv.second);
				
				// Check how many variants in the current range may be handled.
				while (current_lineno < var_lineno)
				{
					// valid_alts() is O(1).
					if (m_alt_checker->valid_alts(current_lineno).empty())
						skipped_lines.emplace_back(current_lineno);
					else
						++current_range_variant_count;
					
					++current_lineno;
				}
				assert(var_lineno == current_lineno);
				
				if (m_min_path_length <= current_range_variant_count)
				{
					// Mark skipped lines in a bit vector.
					std::size_t const total_lines(current_lineno - current_start_lineno + 1);
					sdsl::bit_vector skipped_lines_b(total_lines, 0);
					for (auto const skipped_line : skipped_lines)
					{
						std::size_t const line_idx(skipped_line - current_start_lineno);
						skipped_lines_b[line_idx] = 1;
					}
					
					graph_ranges.emplace_back(
						current_start - buffer_start,
						next_subgraph_start - current_start,
						current_start_lineno,
						current_range_variant_count,
						std::move(skipped_lines_b)
					);
					current_range_variant_count = 0;
					current_start = next_subgraph_start;
					current_start_lineno = current_lineno;
					skipped_lines.clear();
				}
			}
		
			// Add the final graph_range to cover the remaining range.
			auto const lineno_limit(1 + m_record_count + m_reader.last_header_lineno());
			while (current_lineno < lineno_limit)
			{
				if (m_alt_checker->valid_alts(current_lineno).empty())
					skipped_lines.emplace_back(current_lineno);
				else
					++current_range_variant_count;
				
				++current_lineno;
			}
			assert(lineno_limit == current_lineno);
			
			// Mark skipped lines in a bit vector.
			std::size_t const total_lines(current_lineno - current_start_lineno + 1);
			sdsl::bit_vector skipped_lines_b(total_lines, 0);
			for (auto const skipped_line : skipped_lines)
			{
				std::size_t const line_idx(skipped_line - current_start_lineno);
				skipped_lines_b[line_idx] = 1;
			}
			
			graph_ranges.emplace_back(
				current_start - buffer_start,
				buffer_end - current_start,
				current_start_lineno,
				current_range_variant_count,
				std::move(skipped_lines_b)
			);
		}
		
		// Allocate enough space for the subgraphs.
		auto const range_count(graph_ranges.size());
		always_assert(range_count);
		m_subgraphs.resize(range_count);
		m_subgraph_bitmap.resize(range_count);
		m_path_matchings.resize(range_count - 1);
		m_remaining_merge_tasks = range_count - 1;
		
		// Calculate the number of steps.
		m_progress_counter.calculate_step_count(m_alt_checker->records_with_valid_alts(), m_generated_path_count, m_remaining_merge_tasks);
		
		// Update status.
		m_status_logger->log([range_count](){
			std::cerr << "Split the variants into " << range_count << " subgraphs." << std::endl;
		});
		m_status_logger->log_message_progress_bar("Reducing samples…");
		
		// Start each task.
		std::size_t task_idx(0);
		for (auto &range : graph_ranges)
		{
			auto const st(dispatch_semaphore_wait(*m_subgraph_semaphore, DISPATCH_TIME_FOREVER));
			always_assert(0 == st);
			
			auto *task(new read_subgraph_variants_task);
			std::unique_ptr <class task> task_ptr(task);
			init_read_subgraph_variants_task(
				*task,
				task_idx++,
				std::move(range)
			);
			
			m_delegate->store_and_execute(std::move(task_ptr));
		}
	}
}
