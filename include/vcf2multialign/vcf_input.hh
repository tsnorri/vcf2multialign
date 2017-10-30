/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VCF_INPUT_HH
#define VCF2MULTIALIGN_VCF_INPUT_HH

#include <vcf2multialign/file_handling.hh>


namespace vcf2multialign {
	
	class vcf_reader;
	
	
	struct vcf_input
	{
		virtual ~vcf_input() {};
		virtual bool getline(std::string_view &dst) = 0;
		virtual void store_first_variant_offset() = 0;
		virtual void reset_to_first_variant_offset() {}
		virtual void fill_buffer(vcf_reader &reader) = 0;
	};


	class vcf_stream_input_base : public vcf_input
	{
	protected:
		std::string					m_buffer;
		std::istream::pos_type		m_first_variant_offset{0};
		std::size_t					m_len{0};
		std::size_t					m_pos{0};
	
	public:
		vcf_stream_input_base():
			m_buffer(128, '\0')
		{
		}
		
		virtual ~vcf_stream_input_base() = 0;
	
		virtual bool getline(std::string_view &dst) override;
		virtual void store_first_variant_offset() override;
		virtual void reset_to_first_variant_offset() override;
		virtual void fill_buffer(vcf_reader &reader) override;
		
	protected:
		virtual void stream_reset() = 0;
		virtual bool stream_getline() = 0;
		virtual std::size_t stream_tellg() const = 0;
		virtual void stream_read(char *data, std::size_t len) = 0;
		virtual std::size_t stream_gcount() const = 0;
		virtual bool stream_eof() const = 0;
	};
	
	
	template <typename t_stream>
	class vcf_stream_input : public vcf_stream_input_base
	{
	protected:
		t_stream					m_stream;
		
	public:
		using vcf_stream_input_base::vcf_stream_input_base;
		
		t_stream &input_stream() { return m_stream; }
		
	protected:
		virtual void stream_reset() override
		{
			m_stream->clear();
			m_stream->seekg(m_first_variant_offset);
		}
		
		virtual bool stream_getline() override							{ return std::getline(*m_stream, m_buffer); }
		virtual std::size_t stream_tellg() const override				{ return m_stream->tellg(); }
		virtual void stream_read(char *data, std::size_t len) override	{ m_stream->read(data, len); }
		virtual std::size_t stream_gcount() const override				{ return m_stream->gcount(); }
		virtual bool stream_eof() const override						{ return m_stream->eof(); }
	};


	class vcf_mmap_input : public vcf_input
	{
	protected:
		typedef boost::iostreams::stream <
			boost::iostreams::basic_array_source <char>
		>	input_stream;
		
	protected:
		mmap_handle					m_handle;
		std::size_t					m_pos{0};
		std::size_t					m_first_variant_offset{0};
	
	public:
		mmap_handle &handle() { return m_handle; }
		virtual bool getline(std::string_view &dst) override;
		virtual void store_first_variant_offset() override;
		virtual void fill_buffer(vcf_reader &reader) override;
	};
}

#endif
