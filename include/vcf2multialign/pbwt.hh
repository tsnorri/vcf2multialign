/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PBWT_HH
#define VCF2MULTIALIGN_PBWT_HH

#include <algorithm>			// std::iota
#include <libbio/bits.hh>
#include <libbio/int_matrix.hh>
#include <limits>
#include <map>
#include <numeric>
#include <type_traits>			// std::is_unsigned_v
#include <utility>				// std::swap
#include <vector>


namespace vcf2multialign {

	template <typename t_index, typename t_divergence, typename t_count>
	struct pbwt_context
	{
		typedef t_index			index_type;
		typedef t_divergence	divergence_type;
		typedef t_count			count_type;
		constexpr static inline auto const DIVERGENCE_MAX{std::numeric_limits <divergence_type>::max()};
		constexpr static inline auto const COUNT_MAX{std::numeric_limits <count_type>::max()};

		struct divergence_value
		{
			static_assert(std::is_unsigned_v <divergence_type>);

			divergence_type	value{};

			divergence_value() = default;

			/* implicit */ divergence_value(divergence_type const value_):
				value(value_)
			{
			}

			// Place DIVERGENCE_MAX first; needed to get the equivalence class count in find_cut_positions_lambda_min().
			bool operator<(divergence_value const other) const { return 1 + value < 1 + other.value; }
			/* implicit */ operator divergence_type() const { return value; }
		};

		std::vector <index_type>				permutation;
		std::vector <index_type>				prev_permutation;
		std::vector <divergence_value>			divergence;
		std::vector <divergence_value>			prev_divergence;
		std::map <divergence_value, count_type>	divergence_value_counts;

		explicit pbwt_context(count_type const count);
		void update_divergence(libbio::bit_matrix::const_slice_type const slice, divergence_value const kk);
		void swap_vectors();
	};


	template <typename t_index, typename t_divergence, typename t_count>
	pbwt_context <t_index, t_divergence, t_count>::pbwt_context(count_type const count):
		permutation(count),
		divergence(count, DIVERGENCE_MAX)
	{
		if (count)
		{
			divergence[0] = 0;
			divergence_value_counts[0] = 1;
			if (1 < count)
				divergence_value_counts[DIVERGENCE_MAX] = count - 1;
		}
		std::iota(permutation.begin(), permutation.end(), 0);
	}


	template <typename t_index, typename t_divergence, typename t_count>
	void pbwt_context <t_index, t_divergence, t_count>::update_divergence(libbio::bit_matrix::const_slice_type const slice, divergence_value const kk)
	{
		// Mostly following Algorithm 2 in Efficient haplotype matching and storage using the positional
		// Burrows–Wheeler transform (PBWT).

		// Note that the size of the matrix slice may be greater than that of the permutation b.c. the former must be word-aligned.

		// First count the ones.
		t_count zero_idx{};
		auto one_idx([&] -> t_count {
			t_count one_count{};
			for (auto const word : slice.to_span())
				one_count += libbio::bits::count_bits_set(word);
			return prev_permutation.size() - one_count;
		}());

		// Update the sorted order.
		permutation.resize(prev_permutation.size());
		divergence.resize(prev_divergence.size());
		divergence_value pp{kk + 1};
		divergence_value qq{kk + 1};
		for (t_count ii{}; ii < prev_permutation.size(); ++ii)
		{
			auto const val_idx(prev_permutation[ii]);
			auto const pd(prev_divergence[ii]);

			if (pp < pd)
				pp = pd;

			if (qq < pd)
				qq = pd;

			{
				auto const it(divergence_value_counts.find(pd));
				--it->second;
				if (0 == it->second)
					divergence_value_counts.erase(it);
			}

			if (0 == slice[val_idx])
			{
				++divergence_value_counts[pp];
				permutation[zero_idx] = val_idx;
				divergence[zero_idx] = pp;
				++zero_idx;
				pp = 0;
			}
			else
			{
				++divergence_value_counts[qq];
				permutation[one_idx] = val_idx;
				divergence[one_idx] = qq;
				++one_idx;
				qq = 0;
			}
		}
	}


	template <typename t_index, typename t_divergence, typename t_count>
	void pbwt_context <t_index, t_divergence, t_count>::swap_vectors()
	{
		using std::swap;
		swap(permutation, prev_permutation);
		swap(divergence, prev_divergence);
		permutation.clear();
		divergence.clear();
	}
}

#endif
