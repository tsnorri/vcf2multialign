/*
 * Copyright (c) 2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_GRAPH_TO_SEQUENCES_FILE_HANDLING_HH
#define VCF2MULTIALIGN_VARIANT_GRAPH_TO_SEQUENCES_FILE_HANDLING_HH

#include <libbio/buffered_writer/dispatch_io_channel_buffered_writer.hh>
#include <libbio/file_handling.hh>


namespace vcf2multialign {
	
	typedef libbio::dispatch_io_channel_buffered_writer	output_stream_type;
	
	
	inline void open_output_file(char const *path, output_stream_type &of, libbio::writing_open_mode const mode)
	{
		namespace lb = libbio;
		auto const flags(
			O_WRONLY |
			(mode & lb::writing_open_mode::CREATE		? O_CREAT : 0) |	// Create if requested.
			(mode & lb::writing_open_mode::OVERWRITE	? O_TRUNC : O_EXCL)	// Truncate if OVERWRITE given, otherwise require that the file does not exist.
		);
		of = lb::dispatch_io_channel_buffered_writer(path, flags, S_IRUSR | S_IWUSR, 512 * 1024, dispatch_get_main_queue());
	}
	
	
	inline void open_output_file(std::string const &path, output_stream_type &of, libbio::writing_open_mode const mode)
	{
		vcf2multialign::open_output_file(path.data(), of, mode);
	}
	
	
	inline void open_founder_output_file(std::size_t const idx, output_stream_type &of, libbio::writing_open_mode const mode)
	{
		vcf2multialign::open_output_file(std::to_string(1 + idx), of, mode);
	}
}

#endif
