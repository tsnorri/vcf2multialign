/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/preprocess/preprocess_logger.hh>


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	inline void dispatch_async_main(dispatch_block_t block)
	{
		dispatch_async(dispatch_get_main_queue(), block);
	}
}


namespace vcf2multialign {
	
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
