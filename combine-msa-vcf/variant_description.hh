/*
 * Copyright (c) 2020-2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMBINE_MSA_VARIANT_DESCRIPTION_HH
#define VCF2MULTIALIGN_COMBINE_MSA_VARIANT_DESCRIPTION_HH

#include <libbio/assert.hh>
#include <string>
#include <vector>
#include "types.hh"


namespace vcf2multialign {
	
	enum class variant_origin : std::uint8_t
	{
		MSA,
		VC
	};
	
	
	struct variant_description
	{
	protected:
		struct private_tag {};
		
	public:
		// The genotype is stored in a std::vector <bool> here.
		// A pair of integers could be used instead, though, since
		// at the moment we only handle unphased VCF files and thus
		// the order of the values does not matter.
		
		std::string					ref;
		std::string					ref_src;
		std::string					alt;
		std::string					id; // Currently we allow only one ID.
		std::vector <bool>			genotype;
		std::vector <std::string>	filters;
		std::size_t					position{};
		std::int32_t				overlap_count{};
		variant_origin				origin{};
		bool						is_skipped{};
		bool						had_alt_eq_to_ref{};
		
	protected:
		template <typename t_ref_string, typename t_alt_string>
		variant_description(
			std::size_t const position_,
			t_ref_string &&ref_,
			t_alt_string &&alt_,
			std::size_t const ploidy,
			std::int32_t const overlap_count_,
			variant_origin const origin_,
			private_tag const
		):
			ref(std::forward <t_ref_string>(ref_)),
			alt(std::forward <t_alt_string>(alt_)),
			genotype(ploidy, false),
			position(position_),
			overlap_count(overlap_count_),
			origin(origin_)
		{
			libbio_assert_lte(overlap_count, ploidy);
			// Assign the GT values.
			genotype.resize(overlap_count);
			genotype.resize(ploidy, true);
		}
		
	public:
		variant_description() = default;
		
		variant_description(std::size_t const ploidy, std::int32_t const gt_count, variant_origin const origin_):
			genotype(ploidy, true),
			origin(origin_)
		{
			libbio_assert_lte(gt_count, ploidy);
			// Assign the GT values.
			genotype.resize(gt_count);
			genotype.resize(ploidy, false);
		}
		
		template <typename t_ref_string, typename t_alt_string>
		variant_description(
			std::size_t const position_,
			t_ref_string &&ref_,
			t_alt_string &&alt_,
			std::size_t const ploidy,
			std::int32_t const overlap_count_,
			variant_origin const origin_
		):
			variant_description(
				position_,
				std::forward <t_ref_string>(ref_),
				std::forward <t_alt_string>(alt_),
				ploidy,
				(overlap_count_ < 0 ? overlap_count_ : std::min <std::int32_t>(overlap_count_, ploidy)),
				origin_,
				private_tag()
			)
		{
		}
	};


	inline std::ostream &operator<<(std::ostream &os, variant_origin const origin)
	{
		switch (origin)
		{
			case variant_origin::MSA:
				os << "MSA";
				break;

			case variant_origin::VC:
				os << "VC";
				break;
		}

		return os;
	}


	inline std::ostream &operator<<(std::ostream &os, variant_description const &desc)
	{
		os
			<< "position: " << desc.position
			<< " id: " << desc.id
			<< " ref: " << desc.ref
			<< " ref_src: " << desc.ref_src
			<< " alt: " << desc.alt
			<< " origin: " << desc.origin;
		return os;
	}
}

#endif
