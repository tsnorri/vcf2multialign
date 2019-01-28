/*
 * Copyright (c) 2017-2018 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/format.hpp>
#include <sys/stat.h>
#include <libbio/fasta_reader.hh>
#include <libbio/vector_source.hh>
#include <vcf2multialign/read_single_fasta_seq.hh>
#include <vcf2multialign/types.hh>

namespace lb	= libbio;
namespace v2m	= vcf2multialign;


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
			lb::vector_source <v2m::vector_type> &vs
		)
		{
			std::cerr << "read sequence of length " << seq_length << '.' << std::endl;
			
			// Use ADL.
			using std::swap;
			seq->resize(seq_length);
			swap(*reference, *seq);
			vs.put_vector(seq);
			++i;
		}
		
		void start() {}
		void finish() { libbio_always_assert_msg(1 == i, "Expected to have read exactly one sequence."); }
	};
}


namespace vcf2multialign {
	
	// Read the contents of a FASTA file into a single sequence.
	void read_single_fasta_seq(lb::file_istream &ref_fasta_stream, vector_type &reference)
	{
		typedef lb::vector_source <vector_type> vector_source;
		typedef lb::fasta_reader <vector_source, callback> fasta_reader;
		
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
				libbio_fail(msg.c_str());
			}
			
			// Preallocate space for the reference.
			std::cerr << "Preallocating a vector of size " << sb.st_size << "…";
			std::unique_ptr <vector_type> vec_ptr;
			vs.get_vector(vec_ptr);
			vec_ptr->resize(sb.st_size);
			vs.put_vector(vec_ptr);
			std::cerr << " done." << std::endl;
		}
		
		callback cb(reference);
		fasta_reader reader;
		
		std::cerr << "Reading reference FASTA into memory…" << std::endl;
		reader.read_to_vector_from_stream(ref_fasta_stream, vs, cb);
	}
}
