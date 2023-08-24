/*
 * Copyright (c) 2017-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <sys/stat.h>
#include <libbio/fasta_reader.hh>
#include <vcf2multialign/utility/read_single_fasta_seq.hh>
#include <vcf2multialign/types.hh>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace vcf2multialign {
	
	// Read the contents of a FASTA file into a single sequence.
	void read_single_fasta_seq(char const *fasta_path, vector_type &seq, char const *seq_name, bool const logs_status)
	{
		if (logs_status)
			lb::log_time(std::cerr) << "Reading FASTA into memory…";

		auto const res(lb::read_single_fasta_sequence(fasta_path, seq, seq_name));
		if (!res)
		{
			std::cerr << "Unable to read reference sequence";
			if (seq_name)
				std::cerr << " with identifier '" << seq_name << '\'';
			std::cerr << " from FASTA at path '" << fasta_path << ".\n";
			std::exit(EXIT_FAILURE);
		}

		if (logs_status)
			std::cerr << " Done. Sequence length was " << seq.size() << '.' << std::endl;
	}
}
