/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <sys/stat.h>
#include <vcf2multialign/fasta_reader.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/vector_source.hh>

namespace v2m = vcf2multialign;


namespace {
	
	struct callback {
		v2m::vector_type *reference;
		size_t i{0};
		
		callback(v2m::vector_type &reference_):
			reference(&reference_)
		{
		}
		
		void handle_sequence(
			std::string const &identifier,
			std::unique_ptr <v2m::vector_type> &seq,
			size_t const &seq_length,
			v2m::vector_source <v2m::vector_type> &vs
		)
		{
			// Use ADL.
			using std::swap;
			seq->resize(seq_length);
			swap(*reference, *seq);
			vs.put_vector(seq);
		}
		
		void finish() {}
	};
	
	
	class read_context : public v2m::status_logger_delegate
	{
	protected:
		typedef v2m::vector_source <v2m::vector_type>		vector_source;
		typedef v2m::fasta_reader <vector_source, callback>	fasta_reader;
		
	protected:
		fasta_reader m_reader;
		
	public:
		virtual std::size_t step_count() const { return 0; }
		virtual std::size_t current_step() const { return m_reader.current_line(); }
		
		void read_single_fasta_seq(
			v2m::file_istream &ref_fasta_stream,
			v2m::vector_type &reference,
			v2m::status_logger &status_logger
		)
		{
			vector_source vs(1, false);
			
			{
				// Approximate the size of the reference from the size of the FASTA file.
				int const fd(ref_fasta_stream->handle());
				struct stat sb;
				auto const st(fstat(fd, &sb));
				if (0 != st)
				{
					auto const err_str(strerror(errno));
					auto const msg(boost::str(boost::format("Unable to stat the reference file: %s") % err_str));
					v2m::fail(msg.c_str());
				}
				
				// Preallocate space for the reference.
				status_logger.log([size = sb.st_size](){
					std::cerr << "Preallocating a vector of size " << size << "…" << std::flush;
				});
				std::unique_ptr <v2m::vector_type> vec_ptr;
				vs.get_vector(vec_ptr);
				vec_ptr->reserve(sb.st_size);
				vs.put_vector(vec_ptr);
				status_logger.log([](){
					std::cerr << " done." << std::endl;
				}, false);
			}
			
			callback cb(reference);
			status_logger.set_delegate(*this);
			status_logger.log_message_counting("Reading the reference FASTA into memory…");
			m_reader.read_from_stream(ref_fasta_stream, vs, cb);
			status_logger.finish_logging();
			status_logger.log([size = reference.size()](){
				std::cerr << "Read sequence of length " << size << '.' << std::endl;
			});
		}
	};
}


namespace vcf2multialign {
	
	// Read the contents of a FASTA file into a single sequence.
	void read_single_fasta_seq(
		file_istream &ref_fasta_stream,
		vector_type &reference,
		status_logger &status_logger
	)
	{
		read_context ctx;
		ctx.read_single_fasta_seq(ref_fasta_stream, reference, status_logger);
	}
}
