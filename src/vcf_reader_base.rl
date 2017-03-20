/*
 Copyright (c) 2017 Tuukka Norri
 This code is licensed under MIT license (see LICENSE for details).
 */

#include <vcf2multialign/vcf_reader_base.hh>


%% machine vcf_parser;
%% write data;


// We would like to post-process strings that have been read and then pass
// them to the variant object with a simple function call at the call site.
// Solution: Templates can have member function pointers as parameters but of a specific type.
// We convert the pointer to (void (t_record_class::*)(void)) for template instantiation and
// back to t_fnt for calling with t_args.
#define HANDLE_STRING_END(FN, ...)	handle_string_end	<(rc_void_mem_fn)(FN)>(FN, ##__VA_ARGS__)
#define HANDLE_INTEGER_END(FN, ...)	handle_integer_end	<(rc_void_mem_fn)(FN)>(FN, ##__VA_ARGS__)


namespace vcf2multialign {
	
	template <vcf_field t_max_field, typename t_variant>
	class vcf_reader_tpl : public vcf_reader_base
	{
	public:
		using vcf_reader_base::vcf_reader_base;
		
	protected:
		t_variant					m_current_variant;
		
	protected:
		typedef void (t_variant::*rc_void_mem_fn)(void);
		
		template <rc_void_mem_fn t_fn, typename t_fnt, typename ... t_args>
		inline void handle_string_end(t_fnt, t_args ... args);
		
		template <rc_void_mem_fn t_fn, typename t_fnt, typename ... t_args>
		inline void handle_integer_end(t_fnt, t_args ... args);
		
		template <int t_continue, int t_break, typename t_fn>
		inline bool check_max_field(vcf_field const field, int const target, t_fn cb);
		
	public:
		void set_parsed_fields(vcf_field max_field);
		
		void read_header();
		
		template <typename t_fn>
		bool parse(t_fn callback);
	};
	
	
	template <vcf_field t_max_field, typename t_variant>
	template <rc_void_mem_fn t_fn, typename t_fnt, typename ... t_args>
	void vcf_reader_tpl <t_max_field, t_variant>::handle_string_end(t_fnt, t_args ... args)
	{
		std::string_view sv(m_start, m_fsm.p - m_start + 1);
		(m_current_variant->*((t_fnt) t_fn))(sv, std::forward <t_args>(args)...);
	}

	
	template <vcf_field t_max_field, typename t_variant>
	template <rc_void_mem_fn t_fn, typename t_fnt, typename ... t_args>
	void vcf_reader_tpl <t_max_field, t_variant>::handle_integer_end(t_fnt, t_args ... args)
	{
		(m_current_variant->*((t_fnt) t_fn))(m_integer, std::forward <t_args>(args)...);
	}
	
	
	template <vcf_field t_max_field, typename t_variant>
	template <
		int t_continue,
		int t_break,
		typename t_fn
	>
	bool vcf_reader_tpl <t_max_field, t_variant>::check_max_field(vcf_field const field, int const target, t_fn cb)
	{
		if (field <= t_max_field || (t_max_field == vcf_field::VARY && field <= m_max_parsed_field))
			return target;
		
		skip_to_next_nl();
		if (!cb(m_current_variant))
			return t_break;
		
		return t_continue;
	}
	
	
	template <vcf_field t_max_field, typename t_variant>
	void vcf_reader_tpl <t_max_field, t_variant>::set_parsed_fields(vcf_field const max_field)
	{
		if (t_max_field != vcf_field::VARY)
			throw std::runtime_error("t_max_field not set to VARY");
		
		vcf_reader_base::set_parsed_fields(max_field);
	}
	
	
	template <vcf_field t_max_field, typename t_variant>
	void vcf_reader_tpl <t_max_field, t_variant>::read_header()
	{
		vcf_reader_base::read_header(t_max_field);
		
		using std::swap;
		t_variant var(sample_count());
		swap(m_current_variant, var);
	}
	
	
	template <vcf_field t_max_field, typename t_variant>
	template <typename t_fn>
	bool vcf_reader_tpl <t_max_field, t_variant>::parse(t_fn cb)
	{
		typedef t_record_class rc;
		bool retval(true);
		int cs(0);
		// FIXME: add throwing EOF actions?
		%%{
			variable p		m_fsm.p;
			variable pe		m_fsm.pe;
			variable eof	m_fsm.eof;
			
			action start_string {
				m_start = fpc;
			}
			
			action start_integer {
				m_integer = 0;
			}
			
			action update_integer {
				// Shift and add the current number.
				m_integer *= 10;
				m_integer += fc - '0';
			}
			
			action end_sample_field {
				// Check the current field is followed by another field, another sample or a new record.
				switch (fc)
				{
					case ':':
						fgoto sample_rec_f;

					case '\t':
					{
						if (m_format.size() != m_format_idx)
							throw std::runtime_error("Not all fields present in the sample");
						
						m_format_idx = 0;
						fgoto sample_rec_f;
					}
					
					case '\n':
					{
						if (m_format.size() != m_format_idx)
							throw std::runtime_error("Not all fields present in the sample");
						
						cb(m_current_variant);
						
						fgoto main;
					}
					
					default:
						throw std::runtime_error("Unexpected character");
				}
			}
			
			action error {
				throw std::runtime_error("Unexpected character");
			}
			
			tab			= '\t';
			
			id_string	= '<' alnum+ '>';
			
			chrom_id	= (alnum+)
				>(start_string)
				%{ HANDLE_STRING_END(&rc::set_chrom_id); };
			
			pos			= (digit+)
				>(start_integer)
				$(update_integer)
				%{ HANDLE_INTEGER_END(&rc::set_pos); };
			
			id_part		= (alnum+)
				>(start_string)
				%{ HANDLE_STRING_END(&rc::set_id, m_idx++); };
			id_rec		= (id_part (';' id_part)*) >{ m_idx = 0; };
			
			ref			= ([ACGTN]+)
				>(start_string)
				%{ HANDLE_STRING_END(&rc::set_ref); };
			
			# FIXME: add breakends.
			simple_alt	= ([ACGTN]+) %{ m_alt_is_complex = false; };
			complex_alt	= ([ACGTN*]+) %{ m_alt_is_complex = true; };
			alt_part	= (simple_alt | complex_alt)
				>(start_string)
				%{ HANDLE_STRING_END(&rc::set_alt, m_idx++, m_alt_is_complex); };
			alt			= alt_part (',' alt_part)*;
			
			qual		= (digit+)
				>(start_integer)
				$(update_integer)
				%{ HANDLE_INTEGER_END(&rc::set_qual); };
			
			# FIXME: add actions.
			filter_pass	= 'PASS';
			filter_part	= alnum+;
			filter		= (filter_pass | (filter_part (';' filter_part)*));
			
			# FIXME: add actions.
			info_key	= (alnum | '_') +;
			info_str	= (alnum | [_.]) +;
			info_val	= info_str (',' info_str)*;
			info_part	= info_key ('=' info_val)?;
			info		= (info_part (';' info_part)*);
			
			# FIXME: handle other formats.
			format_gt	= ('GT') @{ m_format.emplace_back(format_field::GT); };
			format_dp	= ('DP') @{ m_format.emplace_back(format_field::DP); };
			format_gq	= ('GQ') @{ m_format.emplace_back(format_field::GQ); };
			format_ps	= ('PS') @{ m_format.emplace_back(format_field::PS); };
			format_pq	= ('PQ') @{ m_format.emplace_back(format_field::PQ); };
			format_mq	= ('MQ') @{ m_format.emplace_back(format_field::MQ); };
			format_part	= format_gt | format_dp | format_gq | format_ps | format_pq | format_mq;
			format		= ((format_part (':' format_part)*) >{ m_format.clear(); });
			
			# Genotype field in sample.
			sample_gt_null_allele	= '.'
				%{
					m_integer = NULL_ALLELE;
					HANDLE_INTEGER_END(&rc::set_gt, m_sample_idx - 1, m_idx++, m_gt_is_phased);
				};
				
			sample_gt_allele_idx	= (digit+)
				>(start_integer)
				$(update_integer)
				%{ HANDLE_INTEGER_END(&rc::set_gt, m_sample_idx - 1, m_idx++, m_gt_is_phased); };
				
			sample_gt_part			= (sample_gt_null_allele | sample_gt_allele_idx);
			sample_gt_p				= sample_gt_part ('|' ${ m_gt_is_phased = true; } sample_gt_part)+;
			sample_gt_up			= sample_gt_part ('/' sample_gt_part)+;
			sample_gt				= (sample_gt_part | sample_gt_p | sample_gt_up)
				>{
					m_idx = 0;
					m_gt_is_phased = false;
				};
			
			# FIXME: add actions.
			sample_dp		= (digit+);
			sample_gq		= (digit+);
			sample_ps		= (digit+);
			sample_pq		= (digit+);
			sample_mq		= (digit+);
			
			sep				= '\t';		# Field separator
			ssep			= [\t\n:];	# Sample separator
				
			# Handle a newline and continue.
			main_nl := '\n' %{ fgoto main; }
			break_nl := '\n' % { fbreak; }
			
			# Line start.
			# Apparently main has to be able to read a character, so use fhold.
			# Allow EOF, though, with '?'.
			main := any?
				${
					fhold;
					
					m_idx = 0;
					m_format_idx = 0;
					m_sample_idx = 0;
					++m_lineno;
					m_current_variant.reset();
					
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::CHROM, fentry(chrom_id_f), cb);
				}
				$err(error)
				$eof{ retval = false; };
			
			# #CHROM
			chrom_id_f :=
				(chrom_id sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::POS, fentry(pos_f), cb);
				}
				$err(error);
			
			# POS
			pos_f :=
				(pos sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::ID, fentry(id_f), cb);
				}
				$err(error);
			
			# ID
			id_f :=
				(id_rec sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::REF, fentry(ref_f), cb);
				}
				$err(error);
			
			# REF
			ref_f :=
				(ref sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::ALT, fentry(alt_f), cb);
				}
				$err(error);
			
			# ALT
			alt_f :=
				(alt sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::QUAL, fentry(qual_f), cb);
				}
				$err(error);
			
			# QUAL
			qual_f :=
				(qual sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::FILTER, fentry(filter_f), cb);
				}
				$err(error);
			
			# FILTER
			filter_f :=
				(filter sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::INFO, fentry(info_f), cb);
				}
				$err(error);
			
			# INFO
			info_f :=
				(info sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::FORMAT, fentry(format_f), cb);
				}
				$err(error);
			
			# FORMAT
			# FIXME: try to skip parsing the field and compare to the existing format?
			format_f :=
				(format sep)
				@{
					fgoto *check_max_field <fentry(main_nl), fentry(break_nl)>(vcf_field::ALL, fentry(sample_rec_f), cb);
				}
				$err(error);
			
			# Sample fields
			sample_gt_f := (((sample_gt)) ssep @(end_sample_field)) $err(error);
			sample_dp_f := (((sample_dp)) ssep @(end_sample_field)) $err(error);
			sample_gq_f := (((sample_gq)) ssep @(end_sample_field)) $err(error);
			sample_ps_f := (((sample_ps)) ssep @(end_sample_field)) $err(error);
			sample_pq_f := (((sample_pq)) ssep @(end_sample_field)) $err(error);
			sample_mq_f := (((sample_mq)) ssep @(end_sample_field)) $err(error);
			
			# Sample record
			sample_rec_f := "" >to{
				if (! (m_format_idx < m_format.size()))
					throw std::runtime_error("Format does not match the sample");

				++m_sample_idx;

				// Parse according to the format field.
				auto const idx(m_format_idx);
				++m_format_idx;
				auto const field(m_format[idx]);
				switch (field)
				{
					case format_field::GT:
						fgoto sample_gt_f;
						
					case format_field::DP:
						fgoto sample_dp_f;
						
					case format_field::GQ:
						fgoto sample_gq_f;
						
					case format_field::PS:
						fgoto sample_ps_f;
						
					case format_field::PQ:
						fgoto sample_pq_f;
						
					case format_field::MQ:
						fgoto sample_mq_f;
						
					default:
						throw std::runtime_error("Unexpected format value");
				}
			};
			
			dummy := any
				@{
					// Dummy action,
					// make sure that all the machines are instantiated.
					if (false)
					{
						fgoto chrom_id_f;
						fgoto pos_f;
						fgoto id_f;
						fgoto ref_f;
						fgoto alt_f;
						fgoto qual_f;
						fgoto filter_f;
						fgoto info_f;
						fgoto format_f;
						fgoto sample_rec_f;
						
						fgoto sample_gt_f;
						fgoto sample_dp_f;
						fgoto sample_gq_f;
						fgoto sample_ps_f;
						fgoto sample_pq_f;
						fgoto sample_mq_f;
					}
				};
			
			write init;
			write exec;
		}%%
		
		return retval;
	}
}
