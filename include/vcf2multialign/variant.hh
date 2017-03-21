/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_HH
#define VCF2MULTIALIGN_VARIANT_HH

#include <experimental/string_view>
#include <vcf2multialign/types.hh>
#include <vector>


// XXX Hack.
namespace std {
	using std::experimental::string_view;
}


namespace vcf2multialign {
	
	class vcf_reader;
	
	
	struct genotype
	{
		std::size_t	alt{0};
		bool		is_phased{false};
	};
	
	
	struct sample
	{
		std::vector <genotype>	genotype;
	};
	
	
	class variant_base
	{
		friend class vcf_reader;
		
	protected:
		std::vector <sample>	m_samples;
		std::size_t				m_pos{0};
		std::size_t				m_qual{0};
		std::size_t				m_lineno{0};
		
	public:
		variant_base(std::size_t sample_count):
			m_samples(1 + sample_count)
		{
		}
		
		variant_base(variant_base const &) = default;
		variant_base(variant_base &&) = default;
		variant_base &operator=(variant_base const &) & = default;
		variant_base &operator=(variant_base &&) & = default;
		
		void set_lineno(std::size_t const lineno) { m_lineno = lineno; }
		void set_pos(std::size_t const pos) { m_pos = pos; }
		void set_qual(std::size_t const qual) { m_qual = qual; }
		void set_gt(std::size_t const alt, std::size_t const sample, std::size_t const idx, bool const is_phased);
		void reset() { m_samples.clear(); }	// FIXME: does this cause the genotype vectors to be deallocated?

		size_t lineno() const										{ return m_lineno; }
		size_t pos() const											{ return m_pos; };
		size_t zero_based_pos() const;
		sample const &sample(std::size_t const sample_idx) const	{ return m_samples.at(sample_idx); }
	};
	
	
	template <typename t_string>
	class variant_tpl : public variant_base
	{
		friend class vcf_reader;
		
		template <typename>
		friend class variant_tpl;
		
	protected:
		std::vector <t_string>			m_alts;
		std::vector <t_string>			m_id;
		t_string						m_chrom_id;
		t_string						m_ref;
		
	protected:
		template <typename t_other_string>
		void copy_vectors(variant_tpl <t_other_string> const &other);
		
	public:
		variant_tpl(std::size_t sample_count = 0):
			variant_base(sample_count)
		{
		}

		template <typename t_other_string>
		variant_tpl(variant_tpl <t_other_string> const &other):
			variant_base(other),
			m_chrom_id(other.m_chrom_id),
			m_ref(other.m_ref)
		{
			copy_vectors(other);
		}
		
		template <typename t_other_string>
		bool operator==(variant_tpl <t_other_string> const &other) const;
		
		template <typename t_other_string>
		variant_tpl &operator=(variant_tpl <t_other_string> const &other);
		
		std::vector <t_string> const &alts() const	{ return m_alts; }
		t_string const &ref() const					{ return m_ref; }
		
		void reset() { m_alts.clear(); m_id.clear(); };
		void set_chrom_id(std::string_view const &chrom_id) { m_chrom_id = chrom_id; }
		void set_ref(std::string_view const &ref) { m_ref = ref; }
		void set_id(std::string_view const &id, std::size_t const pos);
		void set_alt(std::string_view const &alt, std::size_t const pos, bool const is_complex);
	};
	
	
	// Transient in the sense that strings point to the working
	// buffer of the reader.
	class transient_variant : public variant_tpl <std::string_view>
	{
		friend class vcf_reader;
		
	public:
		using variant_tpl::variant_tpl;
		void reset();
	};
	
	
	class variant : public variant_tpl <std::string>
	{
		friend class vcf_reader;
		
	public:
		using variant_tpl::variant_tpl;
		void reset();
		
	};
	
	
	template <typename t_string>
	template <typename t_other_string>
	void variant_tpl <t_string>::copy_vectors(variant_tpl <t_other_string> const &other)
	{
		m_alts.reserve(other.m_alts.size());
		m_id.reserve(other.m_id.size());
		
		for (auto const &alt : other.m_alts)
			m_alts.emplace_back(alt);
		
		for (auto const &id : other.m_id)
			m_id.emplace_back(id);
	}
	
	
	template <typename t_string>
	template <typename t_other_string>
	bool variant_tpl <t_string>::operator==(variant_tpl <t_other_string> const &other) const
	{
		return this == &other;
	}
	
	
	template <typename t_string>
	template <typename t_other_string>
	auto variant_tpl <t_string>::operator=(variant_tpl <t_other_string> const &other) -> variant_tpl &
	{
		if (*this != other)
		{
			variant_base::operator=(other);
			m_chrom_id = other.m_chrom_id;
			m_ref = other.m_ref;
			copy_vectors(other);
		}
		return *this;
	}
	
	
	template <typename t_string>
	void variant_tpl <t_string>::set_id(std::string_view const &id, std::size_t const pos)
	{
		if (! (pos < m_id.size()))
			m_id.resize(pos + 1);
		
		m_id[pos] = id;
	}
	
	
	template <typename t_string>
	void variant_tpl <t_string>::set_alt(std::string_view const &alt, std::size_t const pos, bool const is_complex)
	{
		if (is_complex)
			throw std::runtime_error("Only simple ALTs are handled");
		
		if (! (pos < m_alts.size()))
			m_alts.resize(pos + 1);
		
		m_alts[pos] = alt;
	}
}

#endif
