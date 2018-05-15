/*
 Copyright (c) 2017-2018 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <experimental/optional>
#include <vcf2multialign/error_logger.hh>


#ifdef __GNUC__
#   if __GNUC__ < 7
// XXX Hack.
namespace std {
	using std::experimental::optional;
	using std::experimental::nullopt;
}
#	endif
#endif


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	void log(
		std::ostream &output_stream,
		char const *reason,
		std::size_t const line_1,
		std::optional <std::size_t> line_2					= std::nullopt,
		std::optional <std::size_t> ref_pos					= std::nullopt,
		std::optional <std::size_t> alt_idx					= std::nullopt,
		std::optional <lb::sv_type> svt						= std::nullopt,
		std::optional <std::size_t> sample_no				= std::nullopt,
		std::optional <std::size_t> chr_idx					= std::nullopt,
		std::optional <std::size_t> handled_with_alt		= std::nullopt,
		std::optional <std::size_t> total_with_alt			= std::nullopt,
		std::optional <std::size_t> handled_with_any_alt	= std::nullopt,
		std::optional <std::size_t> total_with_any_alt		= std::nullopt
	)
	{
		output_stream << line_1 << '\t';
		
		if (line_2)					output_stream << *line_2;
		output_stream << '\t';
		
		if (ref_pos)				output_stream << *ref_pos;
		output_stream << '\t';
		
		if (alt_idx)				output_stream << *alt_idx;
		output_stream << '\t';
		
		if (svt)					output_stream << to_string(*svt);
		output_stream << '\t';
		
		if (sample_no)				output_stream << *sample_no;
		output_stream << '\t';
		
		if (chr_idx)				output_stream << *chr_idx;
		output_stream << '\t';
		
		if (handled_with_alt)		output_stream << *handled_with_alt;
		output_stream << '\t';
		
		if (total_with_alt)			output_stream << *total_with_alt;
		output_stream << '\t';
		
		if (handled_with_any_alt)	output_stream << *handled_with_any_alt;
		output_stream << '\t';
		
		if (total_with_any_alt)		output_stream << *total_with_any_alt;
		output_stream << '\t';
		
		if (reason)					output_stream << reason;
		
		output_stream << "\n";
	}
}


namespace vcf2multialign {
	
	void error_logger::write_header()
	{
		if (is_logging_errors())
		{
			m_output_stream
			<< "VARIANT_1_LINE\t"
			<< "VARIANT_2_LINE\t"
			<< "REF_POS\t"
			<< "ALT_IDX\t"
			<< "SV_TYPE\t"
			<< "SAMPLE\t"
			<< "CHR\t"
			<< "HANDLED_WITH_ALT\t"
			<< "TOTAL_WITH_ALT\t"
			<< "HANDLED_WITH_ANY_ALT\t"
			<< "TOTAL_WITH_ANY_ALT\t"
			<< "REASON\n";
		}
	}
	
	
	void error_logger::log_no_supported_alts(std::size_t const line)
	{
		if (is_logging_errors())
		{
			char const *reason("No supported ALTs");
			log(m_output_stream, reason, line);
		}
	}
	
	
	void error_logger::log_skipped_structural_variant(std::size_t const line, std::size_t const alt_idx, lb::sv_type const svt)
	{
		if (is_logging_errors())
		{
			char const *reason("Skipped structural variant");
			log(m_output_stream, reason, line, std::nullopt, std::nullopt, alt_idx, svt);
		}
	}
	
	
	void error_logger::log_invalid_alt_seq(std::size_t const line, std::size_t const alt_idx, std::string const &alt)
	{
		if (is_logging_errors())
		{
			char const *reason("Unexpected character in ALT");
			log(m_output_stream, reason, line, std::nullopt, std::nullopt, alt_idx);
		}
	}
	
	
	void error_logger::log_conflicting_variants(std::size_t const line_1, std::size_t const line_2)
	{
		if (is_logging_errors())
		{
			char const *reason("Conflicting variants");
			log(m_output_stream, reason, line_1, line_2);
		}
	}
	
	
	void error_logger::log_ref_mismatch(std::size_t const lineno, std::size_t const diff_pos)
	{
		if (is_logging_errors())
		{
			char const *reason("REF does not match the reference sequence (output anyway using the reference)");
			log(m_output_stream, reason, lineno, std::nullopt, diff_pos);
		}
	}

	
	void error_logger::log_overlapping_alternative(
		std::size_t const lineno,
		std::size_t const sample_no,
		std::size_t const chr_idx,
		sample_count const &alt_counts,
		sample_count const &non_ref_total_counts
	)
	{
		char const *reason("Overlapping alternative");
		log(
			m_output_stream,
			reason,
			lineno,
			std::nullopt,
			std::nullopt,
			std::nullopt,
			std::nullopt,
			sample_no,
			chr_idx,
			alt_counts.handled_count,
			alt_counts.total_count,
			non_ref_total_counts.handled_count,
			non_ref_total_counts.total_count
		);
	}
}
