/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_READ_SINGLE_FASTA_STREAM_HH
#define VCF2MULTIALIGN_READ_SINGLE_FASTA_STREAM_HH

#include <libbio/file_handling.hh>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	void read_single_fasta_seq(libbio::file_istream &ref_fasta_stream, vector_type &reference);
}

#endif
