# vcf2multialign

Create multiply-aligned haplotype sequences from a variant call file and a reference sequence.

## Build/Runtime Requirements

- [git-remote-hg](https://github.com/felipec/git-remote-hg) (to clone Lemon)
- [Lemon](http://lemon.cs.elte.hu/trac/lemon) (provided as a Git submodule)

On Linux the following libraries are required:

- [libdispatch](http://nickhutchinson.me/libdispatch/) (provided as a Git submodule)
- [libpthread_workqueue](https://github.com/mheily/libpwq) (provided as a Git submodule)
- [libkqueue](https://github.com/mheily/libkqueue)

## Build Requirements

- Reasonably new compilers for C and C++, e.g. GCC 6 or Clang 3.9. C++17 support is required.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [Ragel State Machine Compiler](http://www.colm.net/open-source/ragel/) (tested with version 6.7)
- [CMake](http://cmake.org)
- [Boost](http://www.boost.org)

## Building

### Short version

1. `git clone https://github.com/tsnorri/vcf2multialign.git`
2. `cd vcf2multialign`
3. `git submodule update --init --recursive`
4. `cp linux-static.local.mk local.mk`
5. Edit local.mk
6. `make -j4`

### Long version

1. Clone the repository with `git clone https://github.com/tsnorri/vcf2multialign.git`.
2. Change the working directory with `cd vcf2multialign`.
3. Run `git submodule update --init --recursive`. This clones the missing submodules and updates their working tree.
4. Create the file `local.mk`. `linux-static.local.mk` is provided as an example and may be copied with `cp linux-static.local.mk local.mk`
5. Edit `local.mk` in the repository root to override build variables. Useful variables include `CC`, `CXX`, `RAGEL` and `GENGETOPT` for C and C++ compilers, gengetopt and Ragel respectively. `BOOST_INCLUDE` is used as preprocessor flags when Boost is required. `BOOST_LIBS` and `LIBDISPATCH_LIBS` are passed to the linker. See `common.mk` for additional variables.
6. Run make with a suitable numer of parallel jobs, e.g. `make -j4`

Useful make targets include:

<dl>
	<dt>all</dt>
	<dd>Build everything</dd>
	<dt>clean</dt>
	<dd>Remove build products except for dependencies (in the <code>lib</code> folder).</dd>
	<dt>clean-all</dt>
	<dd>Remove all build products.</dd>
</dl>


## Running

The tool takes a Variant Call Format file and a FASTA reference file as its inputs. It then proceeds to read the reference into memory and process the variant file. For each chromosome in the samples part of the VCF, a file is opened in the current working directory and a multiply-aligned haplotype sequence is output. Since the number of files opened may exceed user limits, the VCF is processed in multiple passes.

The FASTA file should contain one sequence only. Currently the VCF parser accepts only a subset of all possible VCF files.

Please see `src/vcf2multialign --help` for command line options.
