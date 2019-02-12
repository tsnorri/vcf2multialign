/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_HANDLER_BASE_HH
#define VCF2MULTIALIGN_VARIANT_HANDLER_BASE_HH

#include <libbio/vcf_reader.hh>
#include <set>
#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class variant_handler_base
	{
	protected:
		std::set <std::uint8_t>							m_valid_alts;
		error_logger									*m_error_logger{};
		sv_handling										m_sv_handling_method{};
		
	protected:
		variant_handler_base() = default;
		
		variant_handler_base(
			sv_handling const sv_handling_method,
			error_logger &error_logger
		):
			m_error_logger(&error_logger),
			m_sv_handling_method(sv_handling_method)
		{
		}
		
		bool check_alt_seq(std::string const &alt) const;
		void fill_valid_alts(libbio::variant const &var);
	};
}

#endif
