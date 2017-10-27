/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_FILE_HANDLING_HH
#define VCF2MULTIALIGN_FILE_HANDLING_HH

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <vcf2multialign/cxx_compat.hh>
#include <vcf2multialign/channel_sink.hh>


namespace vcf2multialign {
	
	// Manage a memory-mapped file.
	// Ragel works with pointer ranges rather than streams, so using Boost's mapped_file_sink would add complexity.
	class mmap_handle
	{
	protected:
		char		*m_content;
		std::size_t	m_mapped_size{0};
		
	public:
		mmap_handle() = default;
		mmap_handle(mmap_handle &&other) = default;
		mmap_handle(mmap_handle const &other) = delete;
		~mmap_handle();
		
		void open(int fd);
		void open(char const *path);
		void close();
		
		char const *data() const { return m_content; }
		std::size_t size() const { return m_mapped_size; }
		
		mmap_handle &operator=(mmap_handle const &other) = delete;
		mmap_handle &operator=(mmap_handle &&other) & = default;
		operator std::string_view() const { return std::string_view(m_content, m_mapped_size); }
	};
	
	
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_source>	file_istream;
	typedef boost::iostreams::stream <boost::iostreams::file_descriptor_sink>	file_ostream;
	
	typedef boost::iostreams::stream <channel_sink>								channel_ostream;
	
	void handle_file_error(char const *fname);
	void open_file_for_reading(char const *fname, file_istream &stream);
	void open_file_for_writing(char const *fname, file_ostream &stream, bool const should_overwrite);
	void open_file_channel_for_writing(
		char const *fname,
		channel_ostream &stream,
		dispatch_ptr <dispatch_semaphore_t> const &write_semaphore,
		bool const should_overwrite
	);
}

#endif
