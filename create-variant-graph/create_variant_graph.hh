/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_CREATE_VARIANT_GRAPH_CREATE_VARIANT_GRAPH_HH
#define VCF2MULTIALIGN_CREATE_VARIANT_GRAPH_CREATE_VARIANT_GRAPH_HH


namespace vcf2multialign {
	
	void create_variant_graph(
		char const *reference_file_path,
		char const *variant_file_path,
		char const *cut_position_file_path,
		char const *output_graph_path,
		char const *reference_seq_name,
		char const *chr_name,
		bool const should_overwrite_files
	);
}

#endif
