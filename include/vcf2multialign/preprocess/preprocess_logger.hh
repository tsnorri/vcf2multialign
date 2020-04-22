/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_PREPROCESS_LOGGER_HH
#define VCF2MULTIALIGN_PREPROCESS_PREPROCESS_LOGGER_HH

#include <atomic>
#include <libbio/dispatch.hh>
#include <libbio/file_handling.hh>
#include <libbio/vcf/variant.hh>
#include <mutex>
#include <set>
#include <string>
#include <vcf2multialign/preprocess/variant_partitioner.hh>



namespace vcf2multialign {
	
	class preprocess_logger : public variant_partitioner_delegate
	{
	protected:
		libbio::dispatch_ptr <dispatch_queue_t>			m_queue;
		libbio::file_ostream							m_log_output_stream;
		std::set <std::pair <std::size_t, std::size_t>>	m_reported_overlaps;
		std::mutex										m_mutex{};
		std::atomic_size_t								m_chrom_id_mismatches{}; // We expect this to be accessed not very frequently from different threads.
		
	public:
		void open_log_file(char const *path, bool const should_overwrite_files);
		
		void variant_processor_found_variant_with_chrom_id_mismatch(libbio::vcf::transient_variant const &var) override { ++m_chrom_id_mismatches; }
		void variant_processor_no_field_for_identifier(std::string const &identifier) override;
		void variant_processor_found_variant_with_position_greater_than_reference_length(libbio::vcf::transient_variant const &var) override;
		void variant_processor_found_variant_with_no_suitable_alts(libbio::vcf::transient_variant const &var) override;
		void variant_processor_found_filtered_variant(libbio::vcf::transient_variant const &var, libbio::vcf::info_field_base const &field) override;
		void variant_processor_found_variant_with_ref_mismatch(libbio::vcf::transient_variant const &var, std::string_view const &ref_sub) override;
		void sample_sorter_found_overlapping_variant(libbio::vcf::variant const &var, std::size_t const sample_idx, std::size_t const prev_end_pos) override;

		std::size_t chrom_id_mismatches() const { return m_chrom_id_mismatches; }
		
	protected:
		template <typename t_extra>
		void log_wt(std::size_t const lineno, std::size_t const sample_idx, char const *reason, t_extra const &extra);
		
		template <typename t_extra>
		void log_wt(std::size_t const lineno, char const *reason, t_extra const &extra);
		
		void log_wt(std::size_t const lineno, std::size_t const sample_idx, char const *reason) { log_wt(lineno, sample_idx, reason, '.'); }
		void log_wt(std::size_t const lineno, char const *reason) { log_wt(lineno, reason, '.'); }
	};
}

#endif
