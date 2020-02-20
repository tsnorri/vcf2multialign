/*
 * Copyright (c) 2020 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <libbio/assert.hh>
#include <libbio/cxxcompat.hh>
#include <libbio/file_handling.hh>
#include <libbio/vcf/variant.hh>
#include <libbio/vcf/variant_printer.hh>
#include <libbio/vcf/vcf_reader.hh>
#include <libbio/vcf/vcf_subfield.hh>
#include <range/v3/all.hpp>
#include <vcf2multialign/utility/forwarder.hh>
#include "cmdline.h"
#include "vcf_record_generator.hh"


namespace lb	= libbio;
namespace rsv	= ranges::view;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {

	// Helper for generating an array of objects with types erased from a tuple.
	template <typename t_base, typename t_tuple, std::size_t ... I>
	constexpr auto make_pointers(t_tuple const &values, std::index_sequence <I...>) -> std::array <t_base const *, std::tuple_size_v <t_tuple>>
	{
		return {&std::get <I>(values)...};
	}

	template <typename t_base, typename t_tuple>
	constexpr auto make_pointers(t_tuple const &values)
	{
		return make_pointers <t_base>(values, std::make_index_sequence <std::tuple_size_v <t_tuple>>{});
	}


	// Comparators for the VCF types.
	struct value_comparator {};

	template <typename t_type>
	struct typed_value_comparator : public value_comparator
	{
		virtual bool lhs_is_missing() const = 0;
		virtual bool compare(t_type const &lhs, t_type const &rhs) const = 0;
	};

	template <typename t_type>
	struct equals_value_comparator final : public typed_value_comparator <t_type>
	{
		bool lhs_is_missing() const override { return false; }
		bool compare(t_type const &lhs, t_type const &rhs) const override { return lhs == rhs; }
	};

	// Wrap the comparator into an object that handles either info or genotype values.
	struct field_comparator_base {};

	struct info_field_comparator : public field_comparator_base
	{
		typedef lb::variant_base				container_type;
		typedef lb::vcf_typed_info_field_base	field_base_type;

		template <bool t_is_vector, lb::vcf_metadata_value_type t_value_type>
		using field_tpl = lb::vcf_typed_info_field <t_is_vector, t_value_type>;

		virtual bool lhs_is_missing() const = 0;

		virtual bool compare(
			container_type const &lhs,
			container_type const &rhs,
			field_base_type const &lhs_access,
			field_base_type const &rhs_access
		) const = 0;
	};

	struct genotype_field_comparator : public field_comparator_base
	{
		typedef lb::variant_sample				container_type;
		typedef lb::vcf_typed_info_field_base	field_base_type;

		template <bool t_is_vector, lb::vcf_metadata_value_type t_value_type>
		using field_tpl = lb::vcf_typed_genotype_field <t_is_vector, t_value_type>;

		virtual bool compare(
			container_type const &lhs,
			container_type const &rhs,
			field_base_type const &lhs_access,
			field_base_type const &rhs_access
		) const = 0;
	};

	template <bool t_is_vector, lb::vcf_metadata_value_type t_value_type, typename t_base>
	struct field_comparator final : public t_base
	{
	public:
		typedef typename t_base::field_base_type								field_base_type;
		typedef typename t_base::template field_tpl <t_is_vector, t_value_type>	field_type;
		typedef typename field_type::container_type								container_type;	// Make sure that the name is actually defined.
		typedef typename field_type::value_type									value_type;		// Underlying value type.
	
		static_assert(std::is_same_v <typename t_base::container_type, container_type>);

	protected:
		typed_value_comparator <value_type> const	*m_comparator{};

	public:
		constexpr field_comparator(typed_value_comparator <value_type> const &comparator):
			m_comparator(&comparator)
		{
		}

		virtual bool lhs_is_missing() const { return m_comparator->lhs_is_missing(); }

		bool compare(
			container_type const &lhs,
			container_type const &rhs,
			field_type const &lhs_access,
			field_type const &rhs_access
		) const { return m_comparator->compare(lhs_access(lhs), rhs_access(rhs)); }

		bool compare(
			container_type const &lhs,
			container_type const &rhs,
			field_base_type const &lhs_access,
			field_base_type const &rhs_access
		) const;
	};

	template <bool t_is_vector, lb::vcf_metadata_value_type t_value_type, typename t_base>
	bool field_comparator <t_is_vector, t_value_type, t_base>::compare(
		container_type const &lhs,
		container_type const &rhs,
		field_base_type const &lhs_access,
		field_base_type const &rhs_access
	) const
	{
		return compare(
			lhs,
			rhs,
			dynamic_cast <field_type const &>(lhs_access),
			dynamic_cast <field_type const &>(lhs_access)
		);
	}

	// Helpers for counting.
	constexpr std::size_t g_scalar_type_count(lb::to_underlying(lb::vcf_metadata_value_type::SCALAR_LIMIT) - lb::to_underlying(lb::vcf_metadata_value_type::FIRST));
	constexpr std::size_t g_vector_type_count(lb::to_underlying(lb::vcf_metadata_value_type::VECTOR_LIMIT) - lb::to_underlying(lb::vcf_metadata_value_type::FIRST));

	// Make the comparators (only equals for now).
	template <bool t_is_vector, std::size_t ... I>
	constexpr auto make_value_comparators(std::index_sequence <I...>)
	{
		constexpr auto make_comparator([](auto const J) constexpr {
			constexpr auto const metadata_value_type(static_cast <lb::vcf_metadata_value_type>(J() + lb::to_underlying(lb::vcf_metadata_value_type::FIRST)));
			typedef lb::vcf_value_type_mapping_t <t_is_vector, metadata_value_type> value_type;
			return equals_value_comparator <value_type>();
		});
		return std::make_tuple(make_comparator(std::integral_constant <std::size_t, I>{})...);
	}

	constexpr static auto g_equals_value_comparators(
		std::tuple_cat(
			make_value_comparators <false>(std::make_index_sequence <g_scalar_type_count>()),
			make_value_comparators <true>(std::make_index_sequence <g_vector_type_count>())
		)
	);

	constexpr std::size_t value_comparator_index(bool is_vector, lb::vcf_metadata_value_type value_type)
	{
		return (is_vector * g_scalar_type_count) + lb::to_underlying(value_type) - lb::to_underlying(lb::vcf_metadata_value_type::FIRST);
	}

	// Make the field comparators.
	template <bool t_is_vector, typename t_base, std::size_t ... I>
	constexpr auto make_field_comparators(std::index_sequence <I...>)
	{
		auto make_comparator([](auto const J) constexpr {
			constexpr auto const value_type(static_cast <lb::vcf_metadata_value_type>(J() + lb::to_underlying(lb::vcf_metadata_value_type::FIRST)));
			constexpr auto const &eq_value_cmp(std::get <value_comparator_index(t_is_vector, value_type)>(g_equals_value_comparators));
			typedef field_comparator <t_is_vector, value_type, t_base> return_value_type;
			return return_value_type(eq_value_cmp);
		});
		return std::make_tuple(make_comparator(std::integral_constant <std::size_t, I>{})...);
	}

	constexpr static auto g_info_field_comparators(
		std::tuple_cat(
			make_field_comparators <false, info_field_comparator>(std::make_index_sequence <g_scalar_type_count>()),
			make_field_comparators <true, info_field_comparator>(std::make_index_sequence <g_vector_type_count>())
		)
	);

	// Erase here, too.
	constexpr static auto g_info_field_comparator_pointers(make_pointers <info_field_comparator>(g_info_field_comparators));

	constexpr std::size_t info_field_comparator_index(bool is_vector, lb::vcf_metadata_value_type value_type)
	{
		// Use the same indexing for now.
		return value_comparator_index(is_vector, value_type);
	}


	struct info_field_checker
	{
		info_field_comparator const			*comparator{};
		lb::vcf_typed_info_field_base const	*lhs_field{};
		lb::vcf_typed_info_field_base const	*rhs_field{};

		info_field_checker(
			info_field_comparator const &comparator_,
			lb::vcf_typed_info_field_base const &lhs_field_,
			lb::vcf_typed_info_field_base const &rhs_field_
		):
			comparator(&comparator_),
			lhs_field(&lhs_field_),
			rhs_field(&rhs_field_)
		{
		}

		bool check(lb::variant const &lhs, lb::variant const &rhs) const
		{
			auto const lhs_is_missing(lhs_field->has_value(lhs));
			auto const rhs_is_missing(rhs_field->has_value(rhs));
			if (lhs_is_missing)
			{
				if (rhs_is_missing)
					return true;
				return comparator->lhs_is_missing();
			}
			return comparator->compare(lhs, rhs, *lhs_field, *rhs_field);
		}
	};


	template <bool t_is_superset_variant>
	struct variant_pack
	{
		lb::variant	variant;
		bool		did_succeed{false};
	};


	class variant_checker final : public lb::vcf_reader_delegate
	{
		friend class forwarder <variant_checker>;

	protected:
		vcf_record_generator				m_subset_gen;
		vcf_record_generator				m_superset_gen;
		std::vector <info_field_checker>	m_info_field_checkers;
		lb::variant							m_subset_var;
		bool								m_found_match{true};	// Initially true to handle superset variants with lower POS than the first subset variant.

	public:
		void prepare(char const *subset_path, char const *superset_path, std::span <char *> const &info_key_names);
		void check_contains();

		void vcf_reader_did_parse_metadata(lb::vcf_reader &reader) override { fix_wp_field(reader); }

	protected:
		void prepare_vcf_record_generator(vcf_record_generator &gen, char const *vcf_path);
		void fix_wp_field(lb::vcf_reader &reader);
		void check_found_match() const;
		void handle(variant_pack <false> &&pack);
		void handle(variant_pack <true> &&pack);
		bool compare_to_superset_variant(lb::variant const &var) const;
	};


	void variant_checker::check_found_match() const
	{
		if (!m_found_match)
			std::cerr << "ERROR: Variant in subset line " << m_subset_var.lineno() << " not found in superset.\n";
	}


	void variant_checker::handle(variant_pack <false> &&pack)
	{
		check_found_match();
		m_subset_var = std::move(pack.variant);
		m_found_match = false;
	}

	
	void variant_checker::handle(variant_pack <true> &&pack) // From superset.
	{
		// FIXME: buffer unmatched and log when receiving the next subset variant?
		if (!m_found_match)
			m_found_match = compare_to_superset_variant(pack.variant);
	}


	bool variant_checker::compare_to_superset_variant(lb::variant const &var) const
	{
		// Compare just POS, REF, ALT and the chosen INFO values.
		// FIXME: CHROM?

		if (m_subset_var.zero_based_pos() != var.zero_based_pos())
			return false;

		if (m_subset_var.ref() != var.ref())
			return false;

		auto const &superset_alts(var.alts());
		for (auto const &subset_alt : m_subset_var.alts())
		{
			if (superset_alts.end() == std::find(superset_alts.begin(), superset_alts.end(), subset_alt))
				return false;
		}

		for (auto const &info_checker : m_info_field_checkers)
		{
			if (!info_checker.check(m_subset_var, var))
				return false;
		}

		return true;
	}


	void variant_checker::fix_wp_field(lb::vcf_reader &reader)
	{
		auto &info_fields(reader.metadata().info());
		auto wp_it(info_fields.find("WP"));
		if (wp_it != info_fields.end())
		{
			auto &field(wp_it->second);
			field.set_number(1);
			field.set_value_type(lb::vcf_metadata_value_type::STRING);
		}
	}


	void variant_checker::prepare_vcf_record_generator(vcf_record_generator &gen, char const *vcf_path)
	{
		gen.open_variants_file(vcf_path);
		gen.vcf_reader().set_delegate(*this);
		gen.prepare();
		gen.vcf_reader().set_parsed_fields(lb::vcf_field::INFO);
	}


	void variant_checker::prepare(char const *subset_path, char const *superset_path, std::span <char *> const &info_key_names)
	{
		prepare_vcf_record_generator(m_subset_gen, subset_path);
		prepare_vcf_record_generator(m_superset_gen, superset_path);

		auto const &subset_info_fields(m_subset_gen.vcf_reader().info_fields());
		auto const &superset_info_fields(m_superset_gen.vcf_reader().info_fields());
		for (auto const *key_name_ptr : info_key_names)
		{
			std::string_view const key_name(key_name_ptr, std::strlen(key_name_ptr));
			auto const lhs_it(subset_info_fields.find(key_name));
			auto const rhs_it(superset_info_fields.find(key_name));

			// Check that there are field objects for the given keys.
			if (subset_info_fields.end() == lhs_it)
			{
				std::cerr << "WARNING: No info field found for key '" << key_name << "' in subset.\n";
				continue;
			}

			if (superset_info_fields.end() == rhs_it)
			{
				std::cerr << "WARNING: No info field found for key '" << key_name << "' in superset.\n";
				continue;
			}

			// Check that the fields have a type that we can handle.
			auto const &lhs_access(*lhs_it->second);
			auto const &rhs_access(*rhs_it->second);
			
			if (!lhs_access.uses_vcf_type_mapping())
			{
				std::cerr << "WARNING: Info field for key '" << key_name << "' uses a custom type in subset.\n";
				continue;
			}

			if (!rhs_access.uses_vcf_type_mapping())
			{
				std::cerr << "WARNING: Info field for key '" << key_name << "' uses a custom type in superset.\n";
				continue;
			}

			// Check that the fields’ types match.
			auto const &lhs_access_(dynamic_cast <lb::vcf_typed_info_field_base const &>(lhs_access));
			auto const &rhs_access_(dynamic_cast <lb::vcf_typed_info_field_base const &>(rhs_access));
			auto const lhs_value_type(lhs_access_.get_value_type());
			auto const lhs_uses_vectors(lhs_access_.value_type_is_vector());

			if (rhs_access_.get_value_type() != lhs_value_type)
			{
				std::cerr << "WARNING: Info field value type does not match for key '" << key_name << "'.\n";
				continue;
			}

			if (rhs_access_.value_type_is_vector() != lhs_uses_vectors)
			{
				std::cerr << "WARNING: Info field value cardinality does not match for key '" << key_name << "'.\n";
				continue;
			}

			// Instantiate a comparator.
			auto const idx(info_field_comparator_index(lhs_uses_vectors, lhs_value_type));
			auto const *cmp(g_info_field_comparator_pointers[idx]);
			m_info_field_checkers.emplace_back(*cmp, lhs_access_, rhs_access_);
		}
	}
	
	
	void variant_checker::check_contains()
	{
		auto lrange(
			rsv::generate([this](){
				variant_pack <false> retval;
				retval.did_succeed = m_subset_gen.next_variant(retval.variant);
				return retval;
			})
			| rsv::take_while([](auto const &pack){ return pack.did_succeed; })
			| rsv::move
		);
		auto rrange(
			rsv::generate([this](){
				variant_pack <true> retval;
				retval.did_succeed = m_superset_gen.next_variant(retval.variant);
				return retval;
			})
			| rsv::take_while([](auto const &pack){ return pack.did_succeed; })
			| rsv::move
		);

		forwarder fwd(*this);
		typedef std::tuple <std::size_t, std::uint8_t> proj_return_type;
		ranges::merge(
			lrange,
			rrange,
			fwd,
			ranges::less(),
			[](auto const &pack){ return proj_return_type(pack.variant.zero_based_pos(), 0); },
			[](auto const &pack){ return proj_return_type(pack.variant.zero_based_pos(), 1); }
		);
	}
}


int main(int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		std::exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);  // Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);                  // We don't require any input from the user.

	v2m::variant_checker checker;
	checker.prepare(args_info.subset_arg, args_info.superset_arg, std::span <char *>(args_info.compare_info_key_arg, args_info.compare_info_key_given));
	checker.check_contains();

	return 0;
}
