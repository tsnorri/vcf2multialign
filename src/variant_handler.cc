/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/combine.hpp>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/util.hh>
#include <vcf2multialign/variant_handler.hh>


namespace vcf2multialign {
	
	void variant_handler::handle_variant(variant &var)
	{
		auto const lineno(var.lineno());
		if (0 != m_skipped_variants->count(lineno))
			return;
		
		// Preprocess the ALT field to check that it can be handled.
		if (m_delegate->valid_alts(lineno).empty())
			return;
		
		std::size_t const var_pos(var.zero_based_pos());
		always_assert(var_pos < m_reference->size(), [this, lineno](){
			std::cerr
			<< "Variant position on line " << lineno
			<< " greater than reference length (" << m_reference->size() << ")."
			<< std::endl;
		});
		
		m_delegate->handle_variant(var);
	}
	
	
	void variant_handler::enumerate_sample_genotypes(
		variant const &var,
		// sample_no, chr_idx, alt_idx, is_phased
		std::function <void(std::size_t, uint8_t, uint8_t, bool)> const &cb
	)
	{
		// Enumerate the samples. The item at index zero (REF_SAMPLE_NUMBER) is always haploid.
		std::size_t sample_no(0);
		for (auto const &sample : var.samples())
		{
			// Handle the genotype.
			uint8_t chr_idx(0);
			for (auto const gt : sample.get_genotype())
			{
				auto const alt_idx(gt.alt);
				auto const is_phased(gt.is_phased);
			
				cb(sample_no, chr_idx, alt_idx, is_phased);
				++chr_idx;
			}
			++sample_no;
		}
	}
	
	
	void variant_handler::finish()
	{
		m_error_logger->flush();
		m_delegate->finish();
	}
	
	
	void variant_handler::process_variants()
	{
		auto &reader(m_variant_buffer.reader());
		reader.reset();
		m_delegate->prepare(reader);
		
		dispatch(&m_variant_buffer).async <&variant_buffer::read_input>(*m_parsing_queue);
	}
}
