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
			std::cerr << "Read sequence of length " << seq_length << std::endl;
			
			// Use ADL.
			using std::swap;
			seq->resize(seq_length);
			swap(*reference, *seq);
			vs.put_vector(seq);
		}
		
		void finish() {}
	};
}


namespace vcf2multialign {
	
	// Read the contents of a FASTA file into a single sequence.
	void read_single_fasta_seq(file_istream &ref_fasta_stream, vector_type &reference)
	{
		typedef vector_source <vector_type> vector_source;
		typedef fasta_reader <vector_source, callback> fasta_reader;
		
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
				throw(std::runtime_error(msg));
			}
			
			// Preallocate space for the reference.
			std::cerr << "Preallocating a vector of size " << sb.st_size << "…";
			std::unique_ptr <vector_type> vec_ptr;
			vs.get_vector(vec_ptr);
			vec_ptr->reserve(sb.st_size);
			vs.put_vector(vec_ptr);
			std::cerr << " done." << std::endl;
		}
		
		callback cb(reference);
		fasta_reader reader;
		
		std::cerr << "Reading reference FASTA into memory…" << std::endl;
		reader.read_from_stream(ref_fasta_stream, vs, cb);
	}
}
