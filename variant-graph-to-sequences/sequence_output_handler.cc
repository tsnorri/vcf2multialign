/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include "sequence_output_handler.hh"


namespace lb	= libbio;
namespace v2m	= vcf2multialign;


namespace {
	
	template <typename t_type>
	struct malloc_deleter
	{
		void operator()(t_type *ptr) const { free(ptr); }
	};
	
	template <typename t_type>
	using malloc_ptr = std::unique_ptr <t_type, malloc_deleter <t_type>>;
}


namespace vcf2multialign {
	
	constexpr inline std::size_t const IO_CHANNEL_WRITER_BUFFER_SIZE{512 * 1024};
	
	
	void subprocess_output_adapter::open_output(char const *name, libbio::writing_open_mode const)
	{
		namespace lb = libbio;
		libbio_assert(m_subprocess);
		libbio_assert(m_output_stream);
		
		{
			// Play it safe. We could just const_cast m_exec_path.data() and name
			// since execvp() would like to have a char *const argv[], but
			// it could then modify the char arrays.
			malloc_ptr <char> exec_path(strdup(m_exec_path));
			malloc_ptr <char> name_(strdup(name));
			
			char * const args[]{exec_path.get(), name_.get(), nullptr};
			*m_subprocess = subprocess_type(args);
		}
			
		auto const fd(m_subprocess->stdin_handle().get());
		*m_output_stream = output_stream_type(
			fd,
			IO_CHANNEL_WRITER_BUFFER_SIZE,
			dispatch_get_main_queue(),
			lb::dispatch_io_channel_flags::NONE // m_subprocess owns the file handle.
		);
	}
	
	
	void file_output_adapter::open_output(char const *path, libbio::writing_open_mode const mode)
	{
		namespace lb = libbio;
		libbio_assert(m_output_stream);
		
		lb::file_handle temp_handle(lb::open_file_for_writing(path, mode));
		*m_output_stream = output_stream_type(
			temp_handle.get(),
			IO_CHANNEL_WRITER_BUFFER_SIZE,
			dispatch_get_main_queue(),
			lb::dispatch_io_channel_flags::HAS_RANDOM_ACCESS | lb::dispatch_io_channel_flags::OWNS_FILE_DESCRIPTOR
		);
		temp_handle.release();
	}
	
	
	void subprocess_sequence_output_handler::prepare_outputs(std::size_t const count)
	{
		m_subprocesses.resize(count);
		m_output_streams.resize(count);
	}
	
	
	void subprocess_sequence_output_handler::close_outputs()
	{
		m_output_streams.clear();
		m_subprocesses.clear();
	}
	
	
	void subprocess_sequence_output_handler::process_outputs(range const rng, std::function <void(std::size_t, output_adapter &)> &cb)
	{
		subprocess_output_adapter adapter(m_exec_path);
		libbio_assert_eq(m_subprocesses.size(), m_output_streams.size());
		rng.check(m_subprocesses);
		auto pos(rng.position());
		auto subprocess_it(m_subprocesses.begin() + pos);
		auto os_it(m_output_streams.begin() + pos);
		auto const subprocess_end(subprocess_it + rng.length());
		while (subprocess_it != subprocess_end)
		{
			adapter.assign(*subprocess_it, *os_it);
			cb(pos, adapter);
			++pos;
			++subprocess_it;
			++os_it;
		}
	}
	
	
	void subprocess_sequence_output_handler::process_last_output(std::function <void(output_adapter &)> &cb)
	{
		subprocess_output_adapter oa(m_exec_path);
		oa.assign(m_subprocesses.back(), m_output_streams.back());
		cb(oa);
	}
	
	
	void file_sequence_output_handler::process_outputs(range const rng, std::function <void(std::size_t, output_adapter &)> &cb)
	{
		file_output_adapter adapter;
		rng.check(m_output_streams);
		auto pos(rng.position());
		auto it(m_output_streams.begin() + pos);
		auto const end(it + rng.length());
		while (it != end)
		{
			adapter.assign(*it);
			cb(pos, adapter);
			++pos;
			++it;
		}
	}
	
	
	void file_sequence_output_handler::process_last_output(std::function <void(output_adapter &)> &cb)
	{
		file_output_adapter oa;
		oa.assign(m_output_streams.back());
		cb(oa);
	}
}
