/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_READ_SINGLE_FASTA_STREAM_HH
#define VCF2MULTIALIGN_UTILITY_READ_SINGLE_FASTA_STREAM_HH

#include <libbio/file_handling.hh>
#include <libbio/mmap_handle.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	void read_single_fasta_seq(char const *fasta_path, vector_type &seq, char const *seq_name = nullptr, bool const logs_status = true);
	void read_single_fasta_seq(libbio::mmap_handle <char> &fasta_handle, vector_type &seq, char const *seq_name = nullptr, bool const logs_status = true);
}

#endif
