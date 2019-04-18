/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstddef>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {

	void generate_haplotypes(
		char const *reference_fname,
		char const *variants_fname,
		char const *ref_seq_name,
		char const *chromosome_name,
		output_type const ot,
		char const *out_fname,
		char const *out_reference_fname,
		char const *report_fname,
		char const *null_allele_seq,
		std::size_t const chunk_size,
		bool const should_overwrite_files,
		bool const should_check_ref
	 );
}
