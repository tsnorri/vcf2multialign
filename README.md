# vcf2multialign

Create a set of reference-guided multiple-aligned haplotype sequences from a variant call file and a reference sequence. Please see [releases](https://github.com/tsnorri/vcf2multialign/releases) for pre-built binaries.

## Build/Runtime Requirements

On Linux, [libbsd](https://libbsd.freedesktop.org/) is required.

## Build Requirements

- Reasonably new compilers for C and C++. We use LLVM 8 for our builds. C++17 support is required.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [Ragel State Machine Compiler](http://www.colm.net/open-source/ragel/) (tested with version 6.7)
- [libdispatch](https://github.com/apple/swift-corelibs-libdispatch) for building on Linux (provided as a Git submodule)
- [CMake](http://cmake.org) for building libdispatch on Linux
- [Ninja](https://ninja-build.org) for building libdispatch on Linux
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

Currently some configuration may need to be done manually by editing `local.mk` as described below.

1. Clone the repository with `git clone https://github.com/tsnorri/vcf2multialign.git`.
2. Change the working directory with `cd vcf2multialign`.
3. Run `git submodule update --init --recursive`. This clones the missing submodules and updates their working tree.
4. Create the file `local.mk`. `linux-static.local.mk` is provided as an example and may be copied with `cp linux-static.local.mk local.mk`
5. Edit `local.mk` in the repository root to override build variables. In addition to GNU Make’s built-in variables, some other ones listed in `common.mk` are used. Useful variables include `CC`, `CXX`, `RAGEL` and `GENGETOPT` for C and C++ compilers, Ragel and gengetopt respectively. `BOOST_INCLUDE` is used as preprocessor flags when Boost is required. `BOOST_LIBS` is passed to the linker.
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

Please see `src/vcf2multialign --help` for a summary of the available command line options.

The tool takes a Variant Call Format file and a FASTA reference file as its inputs. It then proceeds to read the reference into memory and process the variant file. For each chromosome in the samples part of the VCF, a file is opened in the current working directory and a multiply-aligned haplotype sequence is output.

The FASTA file should contain one sequence only. Currently the VCF parser accepts only a subset of all possible VCF files.

### Each sample and chromosome copy require a file descriptor per one pass of the variant file

Since one sequence is generated per each sample and chromosome copy, the tool attempts to open as many files before parsing the variant file. If the number of samples is high, the number of files may exceed shell or operating system limits. In order to process the variant file as few times as possible, the maximum number of file descriptors available to processes started by the shell may be changed with `ulimit -n`. For example, `ulimit -n 8192` sets the maximum number to 8192. Another option is to reduce the number of samples to be handled in one pass with vcf2multialign’s `--chunk-size` command line option.
