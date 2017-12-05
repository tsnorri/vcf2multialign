/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GENERATE_CONFIGURATION_HH
#define VCF2MULTIALIGN_GENERATE_CONFIGURATION_HH

#include <boost/optional.hpp>
#include <string>
#include <vcf2multialign/types.hh>


namespace vcf2multialign {
	
	struct generate_configuration
	{
		boost::optional <std::string>	out_reference_fname;
		std::string						null_allele_seq;
		sv_handling						sv_handling_method;
		std::size_t						chunk_size{0};
		std::size_t						min_path_length{0};
		std::size_t						generated_path_count{0};
		
		bool							should_overwrite_files{false};
		bool							should_reduce_samples{false};
		bool							should_print_subgraph_handling{false};
		
		generate_configuration(
			char const *out_reference_fname_,
			char const *null_allele_seq_,
			sv_handling const sv_handling_method_,
			std::size_t const chunk_size_,
			std::size_t const min_path_length_,
			std::size_t const generated_path_count_,
			bool const should_overwrite_files_,
			bool const should_reduce_samples_,
			bool const print_subgraph_handling_
		):
			null_allele_seq(null_allele_seq_),
			sv_handling_method(sv_handling_method_),
			chunk_size(chunk_size_),
			min_path_length(min_path_length_),
			generated_path_count(generated_path_count_),
			should_overwrite_files(should_overwrite_files_),
			should_reduce_samples(should_reduce_samples_),
			should_print_subgraph_handling(print_subgraph_handling_)
		{
			finish_init(out_reference_fname_);
		}

	protected:
		void finish_init(char const *out_reference_fname_)
		{
			if (out_reference_fname_)
				out_reference_fname.emplace(out_reference_fname_);
		}
	};
}

#endif
