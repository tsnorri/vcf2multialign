/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/utility/find_first_matching_variant.hh>

namespace lb	= libbio;
namespace vcf	= libbio::vcf;
namespace v2m	= vcf2multialign;


namespace {
	class chr_id_variant_validator final : public vcf::variant_validator
	{
	protected:
		std::string_view const m_chr_id;
		
	public:
		chr_id_variant_validator(std::string_view const chr_id):
			m_chr_id(chr_id)
		{
		}
		
		vcf::variant_validation_result validate(vcf::transient_variant const &variant) override
		{
			if (variant.chrom_id() == m_chr_id)
				return vcf::variant_validation_result::PASS;
			
			return vcf::variant_validation_result::SKIP;
		}
	};


	class line_number_variant_validator final : public vcf::variant_validator
	{
	protected:
		std::size_t m_lineno{};
		
	public:
		line_number_variant_validator(std::size_t const lineno):
			m_lineno(lineno)
		{
		}
		
		vcf::variant_validation_result validate(vcf::transient_variant const &variant) override
		{
			auto const var_lineno(variant.lineno());
			if (var_lineno < m_lineno)
				return vcf::variant_validation_result::SKIP;
			else if (var_lineno == m_lineno)
				return vcf::variant_validation_result::PASS;
			
			return vcf::variant_validation_result::STOP;
		}
	};


	template <typename t_validator>
	bool find_using_validator(vcf::reader &reader, t_validator &validator)
	{
		auto &old_validator(reader.variant_validator());
		reader.set_variant_validator(validator);
		
		bool did_find(false);
		reader.parse([&did_find](auto const &var){
			// Reached if a the validation result was PASS (e.g. chromosome ID match).
			did_find = true;
			return false;
		});
		
		reader.set_variant_validator(old_validator);
		return did_find;

	}
}


namespace vcf2multialign {

	bool find_first_matching_variant(vcf::reader &reader, std::string_view const chr_id)
	{
		chr_id_variant_validator validator(chr_id);
		return find_using_validator(reader, validator);
	}


	bool find_variant_by_line_number(libbio::vcf::reader &reader, std::size_t const lineno)
	{
		line_number_variant_validator validator(lineno);
		return find_using_validator(reader, validator);
	}
}
