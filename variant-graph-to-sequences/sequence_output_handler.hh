/*
 * Copyright (c) 2022 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_SEQUENCE_OUTPUT_HANDLER_HH
#define VCF2MULTIALIGN_SEQUENCE_OUTPUT_HANDLER_HH

#include <libbio/buffered_writer/dispatch_io_channel_buffered_writer.hh>
#include <libbio/file_handling.hh>
#include <libbio/subprocess.hh>
#include "range.hh"


namespace vcf2multialign {
	
	typedef libbio::dispatch_io_channel_buffered_writer	output_stream_type;
	
	
	struct output_adapter
	{
		virtual ~output_adapter() {}
		virtual void open_output(char const *name, libbio::writing_open_mode const mode) = 0;
	};
	
	
	class subprocess_output_adapter final : public output_adapter
	{
	public:
		typedef vcf2multialign::output_stream_type							output_stream_type;
		typedef libbio::subprocess <libbio::subprocess_handle_spec::STDIN>	subprocess_type;
		
	protected:
		subprocess_type		*m_subprocess{};
		output_stream_type	*m_output_stream{};
		char const			*m_exec_path{};
		
	public:
		explicit subprocess_output_adapter(std::string const &exec_path):
			m_exec_path(exec_path.data())
		{
		}
		
		void open_output(char const *name, libbio::writing_open_mode const) override;
		void assign(subprocess_type &subprocess, output_stream_type &os) { m_subprocess = &subprocess; m_output_stream = &os; }
	};
	
	
	class file_output_adapter final : public output_adapter
	{
	public:
		typedef vcf2multialign::output_stream_type	output_stream_type;
		
	protected:
		output_stream_type	*m_output_stream{};
		
	public:
		void open_output(char const *path, libbio::writing_open_mode const mode) override;
		void assign(output_stream_type &os) { m_output_stream = &os; }
	};
	
	
	struct sequence_output_handler
	{
		typedef vcf2multialign::output_stream_type	output_stream_type;
		typedef std::vector <output_stream_type>	output_stream_vector;
		
		virtual ~sequence_output_handler() {}
		virtual void prepare_outputs(std::size_t const count) = 0;
		virtual void process_outputs(range const rng, std::function <void(std::size_t, output_adapter &)> &cb) = 0;
		void process_outputs(range const rng, std::function <void(std::size_t, output_adapter &)> &&cb) { process_outputs(rng, cb); }
		virtual void process_last_output(std::function <void(output_adapter &)> &cb) = 0;
		void process_last_output(std::function <void(output_adapter &)> &&cb) { process_last_output(cb); }
		virtual output_stream_vector &output_streams() = 0;
		virtual void close_outputs() = 0;
	};
	
	
	class sequence_output_handler_output_guard
	{
	protected:
		sequence_output_handler	*m_output_handler{};
		
	public:
		sequence_output_handler_output_guard(sequence_output_handler &output_handler):
			m_output_handler(&output_handler)
		{
		}
		
		~sequence_output_handler_output_guard() { m_output_handler->close_outputs(); }
	};
	
	
	class subprocess_sequence_output_handler final : public sequence_output_handler
	{
	public:
		typedef libbio::subprocess <libbio::subprocess_handle_spec::STDIN>	subprocess_type;
		typedef std::vector <subprocess_type>								subprocess_vector;
		
	protected:
		subprocess_vector		m_subprocesses;		// Need to be deallocated after the output streams.
		output_stream_vector	m_output_streams;
		std::string				m_exec_path;
		
	public:
		subprocess_sequence_output_handler(char const *exec_path):
			m_exec_path(exec_path)
		{
		}
		
		output_stream_vector &output_streams() override { return m_output_streams; }
		void prepare_outputs(std::size_t const count) override;
		void close_outputs() override;
		void process_outputs(range const rng, std::function <void(std::size_t, output_adapter &)> &cb) override;
		void process_last_output(std::function <void(output_adapter &)> &cb) override;
	};
	
	
	class file_sequence_output_handler final : public sequence_output_handler
	{
	protected:
		output_stream_vector	m_output_streams;
		
	public:
		output_stream_vector &output_streams() override { return m_output_streams; }
		void prepare_outputs(std::size_t const count) override { m_output_streams.resize(count); }
		void close_outputs() override { m_output_streams.clear(); }
		void process_outputs(range const rng, std::function <void(std::size_t, output_adapter &)> &cb) override;
		void process_last_output(std::function <void(output_adapter &)> &cb) override;
	};
}

#endif
