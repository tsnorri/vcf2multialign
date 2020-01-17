/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <numeric>
#include <vcf2multialign/variant_format.hh>
#include "utility.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	std::pair <std::uint16_t, std::uint16_t> count_set_genotype_values(lb::variant const &var, std::uint16_t const alt_idx)
	{
		typedef std::pair <std::uint16_t, std::uint16_t> return_type;
		
		auto const &samples(var.samples());
		libbio_assert(!samples.empty());
		auto const *gt_field(v2m::get_variant_format(var).gt);
		libbio_assert(gt_field);
		auto const &first_sample(samples.front());
		auto const &gt((*gt_field)(first_sample)); // vector of sample_genotype
		
		// If 0 == alt_idx, count all non-zero values. Otherwise, count the instances of the given value.
		if (0 == alt_idx)
		{
			auto const count(std::accumulate(gt.begin(), gt.end(), std::int32_t(0), [](auto const sum, auto const &sample_gt) -> std::int32_t {
				return (0 == sample_gt.alt ? sum : 1 + sum);
			}));
			return return_type(count, gt.size());
		}
		else
		{
			auto const count(std::accumulate(gt.begin(), gt.end(), std::int32_t(0), [alt_idx](auto const sum, auto const &sample_gt) -> std::int32_t {
				return (alt_idx == sample_gt.alt ? 1 + sum : sum);
			}));
			return return_type(count, gt.size());
		}
	}
}
