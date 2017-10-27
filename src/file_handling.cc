/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <fcntl.h>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <vcf2multialign/channel_sink.hh>
#include <vcf2multialign/file_handling.hh>


namespace ios	= boost::iostreams;


namespace vcf2multialign {
	void handle_file_error(char const *fname)
	{
		char const *errmsg(strerror(errno));
		std::cerr << "Got an error while trying to open '" << fname << "': " << errmsg << std::endl;
		exit(EXIT_FAILURE);
	}


	void open_file_for_reading(char const *fname, file_istream &stream)
	{
		int fd(open(fname, O_RDONLY));
		if (-1 == fd)
			handle_file_error(fname);

		ios::file_descriptor_source source(fd, ios::close_handle);
		stream.open(source);
		stream.exceptions(std::istream::badbit);
	}
	
	
	void open_file_for_writing(char const *fname, file_ostream &stream, bool const should_overwrite)
	{
		int fd(0);
		if (should_overwrite)
			fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		else
			fd = open(fname, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		
		if (-1 == fd)
			handle_file_error(fname);
		
		ios::file_descriptor_sink sink(fd, ios::close_handle);
		stream.open(sink);
	}
	
	
	void open_file_channel_for_writing(
		char const *fname,
		channel_ostream &stream,
		dispatch_ptr <dispatch_semaphore_t> const &write_semaphore,
		bool const should_overwrite
	)
	{
		int fd(0);
		if (should_overwrite)
			fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
		else
			fd = open(fname, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
		
		if (-1 == fd)
			handle_file_error(fname);
		
		channel_sink sink;
		sink.open(fd, write_semaphore);
		stream.open(sink);
	}
	
	
	mmap_handle::~mmap_handle()
	{
		if (m_mapped_size)
		{
			assert(m_content);
			if (-1 == munmap(m_content, m_mapped_size))
				perror("munmap");
		}
	}
	
	
	void mmap_handle::close()
	{
		if (m_mapped_size)
		{
			assert(m_content);
			if (-1 == munmap(m_content, m_mapped_size))
				throw std::runtime_error(strerror(errno));
			
			m_mapped_size = 0;
		}
	}
	
	
	void mmap_handle::open(int fd)
	{
		struct stat sb{};
		if (-1 == fstat(fd, &sb))
			throw std::runtime_error(strerror(errno));
		
		if (!S_ISREG(sb.st_mode))
			throw std::runtime_error("Trying to memory map a non-regular file");
		
		m_content = static_cast <char *>(mmap(0, sb.st_size, PROT_READ, MAP_FILE, fd, 0));
		if (MAP_FAILED == m_content)
			throw std::runtime_error(strerror(errno));
		
		m_mapped_size = sb.st_size;
	}
	
	
	void mmap_handle::open(char const *path)
	{
		int fd(::open(path, O_RDONLY));
		open(fd);
		if (-1 == ::close(fd))
			throw std::runtime_error(strerror(errno));
	}
}
