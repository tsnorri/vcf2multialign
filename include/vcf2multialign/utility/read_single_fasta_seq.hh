/*
 * Copyright (c) 2017-2023 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_READ_SINGLE_FASTA_STREAM_HH
#define VCF2MULTIALIGN_UTILITY_READ_SINGLE_FASTA_STREAM_HH

#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	void read_single_fasta_seq(char const *fasta_path, vector_type &seq, char const *seq_name = nullptr, bool const logs_status = true);
}

#endif
