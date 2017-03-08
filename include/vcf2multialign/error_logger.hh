/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_ERROR_LOGGER_HH
#define VCF2MULTIALIGN_ERROR_LOGGER_HH

#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class error_logger
	{
	protected:
		file_ostream		m_output_stream;
		
	public:
		file_ostream &output_stream() { return m_output_stream; }
		bool is_logging_errors() const { return m_output_stream.is_open(); }
		
		void write_header();
		void log_overlapping_alternative(std::size_t const lineno, std::size_t const sample_no, uint8_t const chd_idx, std::size_t const handled_non_ref_samples);
	};
}

#endif
