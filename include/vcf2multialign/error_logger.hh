/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_ERROR_LOGGER_HH
#define VCF2MULTIALIGN_ERROR_LOGGER_HH

#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/file_handling.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class error_logger
	{
	protected:
		file_ostream					m_output_stream;
		dispatch_ptr <dispatch_queue_t>	m_worker_queue;
		
	public:
		void prepare();
		file_ostream &output_stream()	{ return m_output_stream; }
		void flush()					{ dispatch(this).async <&error_logger::flush_wt>(*m_worker_queue); }
		bool is_logging_errors() const	{ return m_output_stream.is_open(); }
		void write_header()				{ dispatch(this).async <&error_logger::write_header_wt>(*m_worker_queue); }
		
		void log_no_supported_alts(std::size_t const line)
		{
			if (is_logging_errors())
			{
				dispatch_async_fn(*m_worker_queue, [this, line](){
					log_no_supported_alts_wt(line);
				});
			}
		}
		
		void log_skipped_structural_variant(std::size_t const line, std::size_t const alt_idx, sv_type const svt)
		{
			if (is_logging_errors())
			{
				dispatch_async_fn(*m_worker_queue, [this, line, alt_idx, svt](){
					log_skipped_structural_variant_wt(line, alt_idx, svt);
				});
			}
		}
		
		void log_invalid_alt_seq(std::size_t const line, std::size_t const alt_idx, std::string const &alt)
		{
			if (is_logging_errors())
			{
				// Copy alt.
				dispatch_async_fn(*m_worker_queue, [this, line, alt_idx, alt](){
					log_invalid_alt_seq_wt(line, alt_idx, alt);
				});
			}
		}
		
		void log_conflicting_variants(std::size_t const line_1, std::size_t const line_2)
		{
			if (is_logging_errors())
			{
				dispatch_async_fn(*m_worker_queue, [this, line_1, line_2](){
					log_conflicting_variants_wt(line_1, line_2);
				});
			}
		}
		
		void log_ref_mismatch(std::size_t const lineno, std::size_t const diff_pos)
		{
			if (is_logging_errors())
			{
				dispatch_async_fn(*m_worker_queue, [this, lineno, diff_pos](){
					log_ref_mismatch_wt(lineno, diff_pos);
				});
			}
		}
		
		void log_overlapping_alternative(
			std::size_t const lineno,
			std::size_t const sample_no,
			std::size_t const chr_idx,
			sample_count const &alt_counts,
			sample_count const &non_ref_total_counts
		)
		{
			if (is_logging_errors())
			{
				dispatch_async_fn(*m_worker_queue, [this, lineno, sample_no, chr_idx, alt_counts, non_ref_total_counts]{
					log_overlapping_alternative_wt(lineno, sample_no, chr_idx, alt_counts, non_ref_total_counts);
				});
			}
		}
			
	protected:
		void flush_wt() { if (is_logging_errors()) m_output_stream.flush(); }
		void write_header_wt();
		void log_no_supported_alts_wt(std::size_t const line);
		void log_skipped_structural_variant_wt(std::size_t const line, std::size_t const alt_idx, sv_type const svt);
		void log_invalid_alt_seq_wt(std::size_t const line, std::size_t const alt_idx, std::string const &alt);
		void log_conflicting_variants_wt(std::size_t const line_1, std::size_t const line_2);
		void log_ref_mismatch_wt(std::size_t const lineno, std::size_t const diff_pos);
		
		void log_overlapping_alternative_wt(
			std::size_t const lineno,
			std::size_t const sample_no,
			std::size_t const chr_idx,
			sample_count const &alt_counts,
			sample_count const &non_ref_total_counts
		);
	};
}

#endif
