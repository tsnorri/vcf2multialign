/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GRAPH_RANGE_HH
#define VCF2MULTIALIGN_GRAPH_RANGE_HH

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <vcf2multialign/util.hh>


// Represent a range of the DAG constructed from a variant file.
namespace vcf2multialign {
	
	class graph_range
	{
	protected:
		// For overriding move (and copy).
		struct skipped_lines_type
		{
			sdsl::bit_vector indices;
			sdsl::rank_support_v5 <1> index_rank_1_support;
			
			skipped_lines_type() = default;
			
			skipped_lines_type(sdsl::bit_vector &&indices_):
				indices(std::move(indices_)),
				index_rank_1_support(&indices)
			{
				index_rank_1_support.set_vector(&indices);
			}
			
			skipped_lines_type(skipped_lines_type &&other):
				indices(std::move(other.indices)),
				index_rank_1_support(std::move(other.index_rank_1_support))
			{
				index_rank_1_support.set_vector(&indices);
			}
			
			skipped_lines_type &operator=(skipped_lines_type &&other)
			{
				indices = std::move(other.indices);
				index_rank_1_support = std::move(other.index_rank_1_support);
				index_rank_1_support.set_vector(&indices);
				return *this;
			}
		};
		
	protected:
		skipped_lines_type	m_skipped_lines;
		std::size_t			m_range_start_offset{0};
		std::size_t			m_range_length{0};
		std::size_t			m_start_lineno{0};
		std::size_t			m_variant_count{0};
		
	public:
		graph_range() = default;
		
		graph_range(
			std::size_t	range_start_offset,
			std::size_t	range_length,
			std::size_t start_lineno,
			std::size_t	variant_count,
			sdsl::bit_vector &&skipped_lines
		) noexcept:
			m_skipped_lines(std::move(skipped_lines)),
			m_range_start_offset(range_start_offset),
			m_range_length(range_length),
			m_start_lineno(start_lineno),
			m_variant_count(variant_count)
		{
		}
		
		graph_range(graph_range &&other) = default;
		graph_range &operator=(graph_range &&other) = default;
		
		std::size_t seq_position(std::size_t const var_lineno) const;
		bool contains_var_lineno(std::size_t const var_lineno) const;
		std::size_t range_start_offset() const { return m_range_start_offset; }
		std::size_t range_length() const { return m_range_length; }
		std::size_t start_lineno() const { return m_start_lineno; }
		std::size_t variant_count() const { return m_variant_count; }
	};
}

#endif
