/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_COMPRESS_VARIANTS_HH
#define VCF2MULTIALIGN_COMPRESS_VARIANTS_HH

#include <boost/container/vector.hpp>
#include <map>
#include <utility>
#include <vector>
#include <vcf2multialign/error_logger.hh>
#include <vcf2multialign/variant_handler.hh>


namespace vcf2multialign {
	
	typedef std::pair <std::size_t, uint8_t> variant_sequence_id;


	class variant_sequence
	{
	protected:
		variant_sequence_id				m_seq_id{};
		std::size_t						m_start_pos_1{0};
		std::size_t						m_end_pos{0};
		std::map <std::size_t, uint8_t>	m_alt_indices{};	// ALT indices by POS.
	
	public:
		variant_sequence() = default;
	
		variant_sequence(
			std::size_t const sample_no,
			uint8_t const chr_idx
		):
			m_seq_id(sample_no, chr_idx)
		{
		}
	
		variant_sequence(variant_sequence_id const &seq_id):
			m_seq_id(seq_id)
		{
		}
	
		std::size_t start_pos_1() const { return m_start_pos_1; }
		std::size_t start_pos() const { return m_start_pos_1 - 1; }
		std::size_t end_pos() const { return m_end_pos; }
		std::size_t length() const { return 1 + m_end_pos - m_start_pos_1; }
		variant_sequence_id const &seq_id() const { return m_seq_id; }
		std::size_t sample_no() const { return m_seq_id.first; }
		uint8_t chr_idx() const { return m_seq_id.second; }
		
		bool get_alt(std::size_t const zero_based_pos, uint8_t &alt_idx) const
		{
			auto it(m_alt_indices.find(zero_based_pos));
			if (m_alt_indices.cend() == it)
				return false;
			
			alt_idx = it->second;
			return true;
		}
	
		void set_start_pos(std::size_t const zero_based_pos) { m_start_pos_1 = 1 + zero_based_pos; }
	
		bool equal_sequences(variant_sequence const &other) const
		{
			if (m_start_pos_1 != other.m_start_pos_1)
				return false;
		
			if (m_alt_indices.size() != other.m_alt_indices.size())
				return false;
		
			return std::equal(
				m_alt_indices.cbegin(),
				m_alt_indices.cend(),
				other.m_alt_indices.cbegin(),
				other.m_alt_indices.cend()
			);
		}
	
		void reset()
		{
			m_alt_indices.clear();
		}
	
		void add_alt(std::size_t const zero_based_pos, uint8_t const alt_idx)
		{
			m_alt_indices[zero_based_pos] = alt_idx;
		}
	
		bool assign_id(variant_sequence_id const &seq_id)
		{
			if (m_start_pos_1)
				return false;
		
			m_seq_id = seq_id;
			return true;
		}
	};


	typedef boost::container::vector <std::map <std::size_t, variant_sequence>> range_map;
	
	
	void compress_variants(
		vcf_reader &reader,
		error_logger &error_logger,
		variant_set const &skipped_variants,
		std::size_t const padding_amt,
		range_map &compressed_ranges
	);
}

#endif
