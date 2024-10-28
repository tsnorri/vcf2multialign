/*
 * Copyright (c) 2023-2024 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <algorithm>
#include <catch2/catch_all.hpp>
#include <cstddef>
#include <iterator>
#include <libbio/int_matrix/int_matrix.hh>
#include <ostream>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>		// rc::prop
#include <sstream>
#include <tuple>
#include <utility>
#include <vcf2multialign/transpose_matrix.hh>
#include <vector>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {

	struct matrix_size
	{
		std::size_t height{};
		std::size_t width{};
		std::size_t one_count{};

		matrix_size() = default;

		matrix_size(std::size_t const height_, std::size_t const width_):
			height(height_),
			width(width_),
			one_count(4096 * width_ * height_ / 3)
		{
		}

		bool operator<(matrix_size const &other) const { return height * width < other.height * other.width; }
	};

	inline std::ostream &operator<<(std::ostream &os, matrix_size const &ms) { os << ms.height << " × " << ms.width; return os; }


	struct position
	{
		std::size_t row{};
		std::size_t col{};

		auto to_tuple() const { return std::make_tuple(row, col); }
		bool operator<(position const &other) const { return to_tuple() < other.to_tuple(); }
	};


	struct test_case {
		matrix_size		size{};
		lb::bit_matrix	input;
		lb::bit_matrix	expected;
	};


	void output_dimensions(std::ostream &os, lb::bit_matrix const &mat)
	{
		os << mat.number_of_rows() << "×" << mat.number_of_columns();
	}


	std::string comparison(lb::bit_matrix const &input, lb::bit_matrix const &expected, lb::bit_matrix const &actual)
	{
		std::stringstream os;
		auto const &actual_values(actual.values());
		auto const &expected_values(expected.values());

		if (auto const word_size(expected_values.word_size()); word_size == actual_values.word_size())
		{
			for (std::size_t i(0); i < word_size; ++i)
			{
				if (expected_values.word_at(i) != actual_values.word_at(i))
					os << "Mismatch at word " << i << '\n';
			}
		}

		os << std::hex;

		os << "Expected: ";
		std::copy(expected_values.word_begin(), expected_values.word_end(), std::ostream_iterator <lb::bit_matrix::word_type>(os, " "));
		os << '\n';
		os << "Actual:   ";
		std::copy(actual_values.word_begin(), actual_values.word_end(), std::ostream_iterator <lb::bit_matrix::word_type>(os, " "));
		os << '\n';

		os << std::dec;

		os << "Original dimensions: ";
		output_dimensions(os, input);
		os << '\n';

		os << "Expected dimensions: ";
		output_dimensions(os, expected);
		os << '\n';

		os << "Actual dimensions:   ";
		output_dimensions(os, actual);
		os << '\n';

		return os.str();
	}
}


namespace rc {

	template <>
	struct Arbitrary <matrix_size>
	{
		static Gen <matrix_size> arbitrary()
		{
			return gen::shrink(
				gen::scale(0.5, gen::withSize([](int size){ // The test runs quite slowly with sizes above 70 or so.
					if (0 == size)
						return gen::just(matrix_size{});

					return gen::build <matrix_size>(
						gen::set(&matrix_size::height, gen::inRange(1, 1 + size)),
						gen::set(&matrix_size::width, gen::inRange(1, 1 + size))
					);
				})),
				[](matrix_size &&sz)
				{
					return seq::iterate(
						std::move(sz),
						[](matrix_size &&sz){
							if (sz.one_count)
								--sz.one_count;
							return std::move(sz);
						}
					);
				}
			);
		}
	};


	template <>
	struct Arbitrary <test_case>
	{
		static Gen <test_case> arbitrary()
		{
			return gen::mapcat(gen::arbitrary <matrix_size>(), [](auto const &size){
				if (!size.height)
					return gen::just(test_case{});

				auto const bit_height(64 * size.height);
				auto const bit_width(64 * size.width);
				auto const count(bit_width * bit_height);
				return gen::map(
					gen::unique <std::vector <position>>(
						count / 3,
						gen::construct <position>(
							gen::inRange(std::size_t(0), bit_height),
							gen::inRange(std::size_t(0), bit_width)
						)
					),
					[size, bit_width, bit_height](auto const &positions){
						test_case retval{
							size,
							lb::bit_matrix(bit_height, bit_width, 0),
							lb::bit_matrix(bit_width, bit_height, 0)
						};

						for (auto const &pos : positions)
						{
							retval.input(pos.row, pos.col) |= 1;
							retval.expected(pos.col, pos.row) |= 1;
						}

						return retval;
					}
				);
			});
		}
	};
}


SCENARIO("transpose_matrix works with simple input (1×2)", "[transpose_matrix]")
{
	GIVEN("a simple 1x2 matrix")
	{
		lb::bit_matrix input(64, 128);
		lb::bit_matrix expected(128, 64);
		input(1, 68) |= 1;
		expected(68, 1) |= 1;

		WHEN("the matrix is transposed")
		{
			auto const actual(v2m::transpose_matrix(input));
			THEN("the matrix equals the expected")
			{
				INFO(comparison(input, expected, actual));
				REQUIRE(expected == actual);
			}
		}
	}
}


SCENARIO("transpose_matrix works with simple input (2×1)", "[transpose_matrix]")
{
	GIVEN("a simple 1x2 matrix")
	{
		lb::bit_matrix input(128, 64);
		lb::bit_matrix expected(64, 128);
		input(68, 1) |= 1;
		expected(1, 68) |= 1;

		WHEN("the matrix is transposed")
		{
			auto const actual(v2m::transpose_matrix(input));
			THEN("the matrix equals the expected")
			{
				INFO(comparison(input, expected, actual));
				REQUIRE(expected == actual);
			}
		}
	}
}


SCENARIO("transpose_matrix works with simple input (2×2)", "[transpose_matrix]")
{
	GIVEN("a simple 2x2 matrix")
	{
		lb::bit_matrix input(128, 128);
		lb::bit_matrix expected(128, 128);
		input(68, 1) |= 1;
		expected(1, 68) |= 1;

		WHEN("the matrix is transposed")
		{
			auto const actual(v2m::transpose_matrix(input));
			THEN("the matrix equals the expected")
			{
				INFO(comparison(input, expected, actual));
				REQUIRE(expected == actual);
			}
		}
	}
}


TEST_CASE(
	"transpose_matrix with arbitrary input",
	"[transpose_matrix]"
)
{
	matrix_size max_size{};

	rc::prop(
		"transpose_matrix works with arbitrary input",
		[&max_size](test_case const &tc){
			RC_TAG(tc.input.number_of_rows(), tc.input.number_of_columns());

			auto const actual(v2m::transpose_matrix(tc.input));
			auto const &actual_values(actual.values());
			auto const &expected_values(tc.expected.values());
			if (actual_values != expected_values)
				RC_FAIL(comparison(tc.input, tc.expected, actual));

			max_size = std::max(max_size, tc.size);

			return true;
		}
	);

	RC_LOG() << "Max. size: " << max_size << " words\n";
}
