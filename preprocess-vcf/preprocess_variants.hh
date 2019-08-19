/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_PREPROCESS_PREPROCESS_VARIANTS_HH
#define VCF2MULTIALIGN_PREPROCESS_PREPROCESS_VARIANTS_HH


namespace vcf2multialign {
	
	void preprocess_variants(
		char const *reference_path,
		char const *variants_path,
		char const *output_variants_path,
		char const *reference_seq_name,
		char const *chr_name,
		std::vector <std::string> const &field_names_for_filter_if_set,
		bool const should_overwrite_files
	);
}

#endif
