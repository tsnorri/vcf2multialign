/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_GZIP_SINK_HH
#define VCF2MULTIALIGN_GZIP_SINK_HH

#include <boost/container/pmr/synchronized_pool_resource.hpp>
#include <boost/iostreams/categories.hpp>
#include <iosfwd>
#include <vcf2multialign/cxx_compat.hh>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/polymorphic_sink.hh>
#include <vcf2multialign/util.hh>
#include <zlib.h>


namespace vcf2multialign {
	
	class gzip_sink_impl final : public polymorphic_sink_impl
	{
	protected:
		typedef boost::container::pmr::synchronized_pool_resource	pool_resource_type;
		typedef std::shared_ptr <pool_resource_type>				pool_resource_ptr;
		
	protected:
		static pool_resource_type			*s_pool_resource;
		
	protected:
		z_stream							m_compression_stream{};
		dispatch_ptr <dispatch_queue_t>		m_compression_queue{};
		dispatch_ptr <dispatch_queue_t>		m_channel_queue{};
		dispatch_ptr <dispatch_io_t>		m_channel{};
		dispatch_ptr <dispatch_semaphore_t>	m_compression_semaphore{};
		dispatch_ptr <dispatch_semaphore_t>	m_writing_semaphore{};
		dispatch_ptr <dispatch_group_t>		m_writing_group{};
		dispatch_ptr <dispatch_group_t>		m_closing_group{};
		std::size_t							m_write_size{};
		std::atomic_bool					m_close_flag{false};
		
	public:
		static void init() { s_pool_resource = new pool_resource_type; }
		void close() override;
		//void set_pool_resource(pool_resource_ptr const &resource) { m_d.pool_resource = resource; }
		void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) override { m_closing_group = group; }
		void open(
			int fd,
			dispatch_ptr <dispatch_semaphore_t> const &compression_semaphore,
			dispatch_ptr <dispatch_semaphore_t> const &writing_semaphore
		);
		std::streamsize write(char const *bytes, std::streamsize const n) override;
		
	protected:
		unsigned char *allocate_input_buffer(std::size_t const n);
		void deallocate_input_buffer(unsigned char *buffer, std::size_t const n);
		unsigned char *allocate_compression_buffer();
		void deallocate_compression_buffer(unsigned char *buffer);
		
		void allocate_and_place_compression_buffer();
		void write_compressed_data_and_replace_buffer(bool get_new_buffer);
	};
	
	
	class gzip_sink final
	{
	public:
		typedef char						char_type;
		struct category :
			boost::iostreams::sink_tag,
			boost::iostreams::closable_tag
		{
		};
		
	protected:
		// FIXME: copying should work with our current implementation but come up with a better way.
		gzip_sink_impl	*m_impl{};
		
	public:
		gzip_sink() = default;
		
		void close() { if (m_impl) m_impl->close(); }
		void set_closing_group(dispatch_ptr <dispatch_group_t> const &group) { m_impl->set_closing_group(group); }
		std::streamsize write(char const *bytes, std::streamsize const n) { return m_impl->write(bytes, n); }
		
		void open(
			int fd,
			dispatch_ptr <dispatch_semaphore_t> const &compression_semaphore,
			dispatch_ptr <dispatch_semaphore_t> const &writing_semaphore
		)
		{
			m_impl = new gzip_sink_impl; m_impl->open(fd, compression_semaphore, writing_semaphore);
		}
	};
}

#endif
