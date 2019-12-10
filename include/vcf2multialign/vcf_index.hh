/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCF_INDEX_HH
#define VCF2MULTIALIGN_VCF_INDEX_HH

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	class vcf_index
	{
	public:
		struct position
		{
			std::size_t	pos{};
			std::size_t	file_offset{};
			
			position() = default;
			position(std::size_t const pos_, std::size_t const file_offset_):
				pos(pos_),
				file_offset(file_offset_)
			{
			}
			
			// For Cereal.
			template <typename t_archive>
			void serialize(t_archive &archive, std::uint32_t const version);
		};
	
		struct chromosome_entry
		{
			std::string				chromosome_identifier;
			std::vector <position>	positions;
			
			chromosome_entry() = default;
			
			chromosome_entry(std::string const &chr_id_):
				chromosome_identifier(chr_id_)
			{
			}
			
			chromosome_entry(std::string &&chr_id_):
				chromosome_identifier(std::move(chr_id_))
			{
			}
			
			position &add_position(std::size_t const pos, std::size_t const file_offset) { return positions.emplace_back(pos, file_offset); }
			
			// For Cereal.
			template <typename t_archive>
			void serialize(t_archive &archive, std::uint32_t const version);
		};
		
	protected:
		std::vector <chromosome_entry>	m_entries;
		
	public:
		void clear() { m_entries.clear(); }
		chromosome_entry &add_entry(std::string const &chrom) { return m_entries.emplace_back(chrom); }
		chromosome_entry &add_entry(std::string_view const &chrom) { return m_entries.emplace_back(std::string(chrom)); }
		chromosome_entry &add_entry() { return m_entries.emplace_back(); }
		
		// For Cereal.
		template <typename t_archive>
		void serialize(t_archive &archive, std::uint32_t const version);
	};
	
	
	template <typename t_archive>
	void vcf_index::serialize(t_archive &archive, std::uint32_t const version)
	{
		// Ignore the version for now.
		archive(m_entries);
	}
	
	
	template <typename t_archive>
	void vcf_index::chromosome_entry::serialize(t_archive &archive, std::uint32_t const version)
	{
		// Ignore the version for now.
		archive(chromosome_identifier);
		archive(positions);
	}
	
	
	template <typename t_archive>
	void vcf_index::position::serialize(t_archive &archive, std::uint32_t const version)
	{
		// Ignore the version for now.
		archive(pos);
		archive(file_offset);
	}
}

#endif
