/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <sys/stat.h>
#include <libbio/fasta_reader.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/types.hh>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	class delegate : public libbio::fasta_reader_delegate
	{
	protected:
		v2m::vector_type *m_reference{};
		char const *m_identifier{};
		bool m_found_seq{};
		
	public:
		delegate(v2m::vector_type &reference, char const *identifier):
			m_reference(&reference),
			m_identifier(identifier)
		{
		}
		
		virtual bool handle_comment_line(lb::fasta_reader &reader, std::string_view const &sv) override
		{
			return true;
		}
		
		virtual bool handle_identifier(lb::fasta_reader &reader, std::string_view const &sv) override
		{
			if (m_found_seq)
				return false;
			
			if ((!m_identifier) || sv == m_identifier)
			{
				m_found_seq = true;
				std::cerr << " reading sequence with identifier " << sv << "…" << std::flush;
			}
			
			return true;
		}
		
		virtual bool handle_sequence_line(lb::fasta_reader &reader, std::string_view const &sv) override
		{
			if (m_identifier && !m_found_seq)
				return false;
			
			std::copy(sv.begin(), sv.end(), std::back_inserter(*m_reference));
			return true;
		}
	};
}


namespace vcf2multialign {
	
	// Read the contents of a FASTA file into a single sequence.
	void read_single_fasta_seq(char const *fasta_path, vector_type &seq, char const *seq_name)
	{
		lb::mmap_handle <char> handle;
		handle.open(fasta_path);
		read_single_fasta_seq(handle, seq, seq_name);
	}
	
	
	// FIXME: log conditionally, add a flag to parameters (before seq_name).
	void read_single_fasta_seq(libbio::mmap_handle <char> &fasta_handle, vector_type &seq, char const *seq_name)
	{
		lb::fasta_reader reader;
		delegate cb(seq, seq_name);
		
		lb::log_time(std::cerr);
		std::cerr << "Reading FASTA into memory…";
		reader.parse(fasta_handle, cb);
		std::cerr << " Done. Sequence length was " << seq.size() << '.' << std::endl;
	}
}
