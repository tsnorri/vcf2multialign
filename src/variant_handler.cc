/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/combine.hpp>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/util.hh>
#include <vcf2multialign/variant_handler.hh>


namespace vcf2multialign {
	
	bool variant_handler::check_alt_seq(std::string const &alt) const
	{
		for (auto const c : alt)
		{
			if (! ('A' == c || 'C' == c || 'G' == c || 'T' == c || 'N' == c))
				return false;
		}
		
		return true;
	}
	
	
	void variant_handler::fill_valid_alts(variant const &var)
	{
		// Check that the alt sequence is something that can be handled.
		m_valid_alts.clear();
		
		// Insert everything if checking has already been done.
		if (!m_check_alts)
		{
			for (std::size_t i(0), count(var.alts().size()); i < count; ++i)
				m_valid_alts.emplace_hint(m_valid_alts.cend(), i);
				
			return;
		}
		
		std::size_t i(0);
		auto const lineno(var.lineno());
		
		if (sv_handling::DISCARD == m_sv_handling_method)
		{
			for (auto const &ref : boost::combine(var.alts(), var.alt_sv_types()))
			{
				++i;
				
				auto const alt_svt(ref.get <1>());
				if (sv_type::NONE != alt_svt)
					continue;
				
				auto const &alt(ref.get <0>());
				if (!check_alt_seq(alt))
				{
					m_error_logger->log_invalid_alt_seq(lineno, i, alt);
					continue;
				}
				
				m_valid_alts.emplace(i);
			}
		}
		else
		{
			for (auto const &ref : boost::combine(var.alts(), var.alt_sv_types()))
			{
				++i;
				
				auto const alt_svt(ref.get <1>());
				switch (alt_svt)
				{
					case sv_type::NONE:
					{
						auto const &alt(ref.get <0>());
						if (check_alt_seq(alt))
							m_valid_alts.emplace(i);
						else
							m_error_logger->log_invalid_alt_seq(lineno, i, alt);
						
						break;
					}
					
					case sv_type::DEL:
					case sv_type::DEL_ME:
						m_valid_alts.emplace(i);
						break;
					
					case sv_type::INS:
					case sv_type::DUP:
					case sv_type::INV:
					case sv_type::CNV:
					case sv_type::DUP_TANDEM:
					case sv_type::INS_ME:
					case sv_type::UNKNOWN:
						m_error_logger->log_skipped_structural_variant(lineno, i, alt_svt);
						break;
					
					default:
						fail("Unexpected structural variant type.");
						break;
				}
			}
		}
	}
	
	
	void variant_handler::handle_variant(variant &var)
	{
		auto const lineno(var.lineno());
		if (0 != m_skipped_variants->count(lineno))
			return;
		
		// Preprocess the ALT field to check that it can be handled.
		this->fill_valid_alts(var);
		if (m_valid_alts.empty())
			return;
		
		size_t const var_pos(var.zero_based_pos());
		always_assert(var_pos < m_reference->size(), [this, lineno](){
			std::cerr
			<< "Variant position on line " << lineno
			<< " greater than reference length (" << m_reference->size() << ")."
			<< std::endl;
		});
		
		m_delegate->handle_variant(var);
	}
	
	
	void variant_handler::enumerate_genotype(
		variant &var,
		std::size_t const sample_no,
		std::function <void(uint8_t, std::size_t, bool)> const &cb
	)
	{
		// Get the sample.
		auto const sample(var.sample(sample_no));
		
		// Handle the genotype.
		uint8_t chr_idx(0);
		for (auto const gt : sample.get_genotype())
		{
			auto const alt_idx(gt.alt);
			auto const is_phased(gt.is_phased);
			
			cb(chr_idx, alt_idx, is_phased);
			
			++chr_idx;
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
