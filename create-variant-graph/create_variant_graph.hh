/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CREATE_VARIANT_GRAPH_CREATE_VARIANT_GRAPH_HH
#define VCF2MULTIALIGN_CREATE_VARIANT_GRAPH_CREATE_VARIANT_GRAPH_HH


namespace vcf2multialign {
	
	void create_variant_graph_preprocessed(
		char const *reference_file_path,
		char const *variant_file_path,
		char const *preprocessing_result_file_path,
		char const *output_graph_path,
		char const *reference_seq_name,
		bool const should_overwrite_files
	);
	
	void create_variant_graph_single_pass(
		char const *reference_file_path,
		char const *variant_file_path,
		char const *output_graph_path,
		char const *log_path,
		char const *reference_seq_name,
		char const *chr_id,
		std::vector <std::string> &&field_names_for_filter_if_set,
		bool const should_overwrite_files
	);
}

#endif
