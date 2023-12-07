/*
 * Copyright (c) 2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <type_traits>
#include <vcf2multialign/transpose_matrix.hh>

namespace lb	= libbio;


namespace vcf2multialign {
	
	constexpr std::uint64_t transpose8x8(std::uint64_t const word)
	{
		// Partially from https://stackoverflow.com/a/41046873/856976
		return (
			  (word & 0x0100'0000'0000'0000) >> 49
			| (word & 0x0201'0000'0000'0000) >> 42
			| (word & 0x0402'0100'0000'0000) >> 35
			| (word & 0x0804'0201'0000'0000) >> 28
			| (word & 0x1008'0402'0100'0000) >> 21
			| (word & 0x2010'0804'0201'0000) >> 14
			| (word & 0x4020'1008'0402'0100) >> 7
			| (word & 0x8040'2010'0804'0201)
			| (word & 0x0080'4020'1008'0402) << 7
			| (word & 0x0000'8040'2010'0804) << 14
			| (word & 0x0000'0080'4020'1008) << 21
			| (word & 0x0000'0000'8040'2010) << 28
			| (word & 0x0000'0000'0080'4020) << 35
			| (word & 0x0000'0000'0000'8040) << 42
			| (word & 0x0000'0000'0000'0080) << 49
		);
	}
	
	
	lb::bit_matrix transpose_matrix(lb::bit_matrix const &mat)
	{
		static_assert(std::is_same_v <std::uint64_t, lb::bit_matrix::value_type>);
		
		// Doing this in place would be quite difficult (esp. for non-rectangular matrices).
		auto const src_nrow(mat.number_of_rows());
		auto const src_ncol(mat.number_of_columns());
		if (0 == src_ncol)
			return lb::bit_matrix{};

		lb::bit_matrix dst(src_ncol, src_nrow);
		
		libbio_assert_eq(0, src_nrow % 64);
		libbio_assert_eq(0, src_ncol % 64);
		auto const src_col_groups(src_ncol / 64);
		auto const src_col_words(src_nrow / 64);
		auto const dst_col_words(src_ncol / 64);
		
		auto const &src_values(mat.values());
		auto &dst_values(dst.values());
		libbio_assert_eq(src_values.size(), dst_values.size());
		
		for (std::size_t src_row_word_idx(0); src_row_word_idx < src_col_words; ++src_row_word_idx)
		{
			for (std::uint8_t src_row_byte_idx(0); src_row_byte_idx < 8; ++src_row_byte_idx)
			{
				// Process 64 columns at a time so that we can fill the destination columns one 64-bit word at a time.
				for (std::size_t src_col_group(0); src_col_group < src_col_groups; ++src_col_group)
				{
					std::uint64_t src_blocks[8]{};
					std::uint64_t dst_blocks[8]{};

					// Pack the src_row_byte_idx-th byte from each column in the group to the 8×8 matrix below.
					for (std::uint8_t block_idx(0); block_idx < 8; ++block_idx)
					{
						for (std::uint8_t src_col_idx_add(0); src_col_idx_add < 8; ++src_col_idx_add)
						{
							auto const src_col_idx(64 * src_col_group + 8 * block_idx + src_col_idx_add);
							auto const src_word_idx(src_col_idx * src_col_words + src_row_word_idx);
							std::uint64_t src_word(src_values.word_at(src_word_idx));
							src_word >>= 8 * src_row_byte_idx;
							src_word &= 0xff;
							src_word <<= 8 * src_col_idx_add;
							src_blocks[block_idx] |= src_word;
						}
						
						dst_blocks[block_idx] = transpose8x8(src_blocks[block_idx]);
					}
					
					for (std::uint8_t block_col_idx(0); block_col_idx < 8; ++block_col_idx)
					{
						// Add src_col_group to get the correct word in the current column.
						auto const dst_idx((64 * src_row_word_idx + 8 * src_row_byte_idx + block_col_idx) * dst_col_words + src_col_group);
						auto &dst_word(dst_values.word_at(dst_idx));
						for (std::uint8_t block_idx(0); block_idx < 8; ++block_idx)
						{
							std::uint64_t word(dst_blocks[block_idx]);
							word >>= 8 * block_col_idx;
							word &= 0xff;
							word <<= 8 * block_idx;
							dst_word |= word;
						}
					}
				}
			}
		}
		
		return dst;
	}
}
