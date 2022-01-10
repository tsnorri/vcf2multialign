# vcf2multialign

Create a set of reference-guided multiple-aligned haplotypes or founder sequences from a variant call file and a reference sequence. Please see [releases](https://github.com/tsnorri/vcf2multialign/releases) for pre-built binaries.

## Academic Use

If you use the software in an academic setting we kindly ask you to cite the following paper.

```TeX
@article{Norri2021FounderReconstruction,
  title = {Founder reconstruction enables scalable and seamless pangenomic analysis},
  keywords = {READ ALIGNMENT, GENOMES, GRAPHS, SET, 113 Computer and information sciences},
  author = {Tuukka Norri and Bastien Cazaux and Saska D{\"o}nges and Daniel Valenzuela and Veli M{\"a}kinen},
  year = {2021},
  month = dec,
  day = {15},
  doi = {10.1093/bioinformatics/btab516},
  url = {https://doi.org/10.1093/bioinformatics/btab516},
  language = {English},
  volume = {37},
  pages = {4611--4619},
  journal = {Bioinformatics},
  issn = {1367-4803},
  publisher = {Oxford University Press},
  number = {24}
}
```

## Build/Runtime Requirements

On Linux, [libbsd](https://libbsd.freedesktop.org/) is required.

## Build Requirements

- Reasonably new compilers for C and C++. We use LLVM 8 for our builds. C++17 support is required.
- [GNU gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [Ragel State Machine Compiler](http://www.colm.net/open-source/ragel/) (tested with version 6.8)
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
6. `make -j16`

### Long version

Currently some configuration may need to be done manually by editing `local.mk` as described below.

1. Clone the repository with `git clone https://github.com/tsnorri/vcf2multialign.git`.
2. Change the working directory with `cd vcf2multialign`.
3. Run `git submodule update --init --recursive`. This clones the missing submodules and updates their working tree.
4. Create the file `local.mk`. `linux-static.local.mk` is provided as an example and may be copied with `cp linux-static.local.mk local.mk`
5. Edit `local.mk` in the repository root to override build variables. In addition to GNU Make’s built-in variables, some other ones listed in `common.mk` are used. Useful variables include `CC`, `CXX`, `RAGEL` and `GENGETOPT` for C and C++ compilers, Ragel and gengetopt respectively. `BOOST_INCLUDE` is used as preprocessor flags when Boost is required. `BOOST_LIBS` is passed to the linker.
6. Run make with a suitable numer of parallel jobs, e.g. `make -j16`

Useful make targets include:

<dl>
	<dt>all</dt>
	<dd>Build everything</dd>
	<dt>clean</dt>
	<dd>Remove build products except for dependencies (in the <code>lib</code> folder).</dd>
	<dt>clean-all</dt>
	<dd>Remove all build products.</dd>
</dl>


## Creating sample or founder sequences

Please use the `--help` option with each of the tools for a summary of the available command line options, as well as examples.

### Creating sample sequences

1. Run `create-variant-graph/create_variant_graph --cut-by-overlap-start`. The tool takes a Variant Call Format file and a FASTA file as its inputs. It outputs the variants as a directed acyclic graph in a binary format.

2. Run `variant-graph-to-sequences/variant_graph_to_sequences --output-samples`. The tool takes a FASTA file and a variant graph generated by `create_variant_graph` as it inputs and outputs either haplotypes or founder sequences, one sequence per file. The files will be written to the current working directory.

### Creating founder sequences

1. Run `preprocess-vcf/preprocess_vcf`. The tool takes a Variant Call Format file and a FASTA file as inputs. It outputs a list of optimal cut positions in binary format in order to minimize the number of paths in each segment of the graph. This information is used when outputting founder sequences.

2. Run `create-variant-graph/create_variant_graph --cut-by-precalculated`. The tool takes a Variant Call Format file, a FASTA file and the position file generated by `preprocess_vcf` as its inputs. It outputs the variants as a directed acyclic graph in a binary format.

3. Run `variant-graph-to-sequences/variant_graph_to_sequences --output-founders`. The tool takes a FASTA file and a variant graph generated by `create_variant_graph` as it inputs and outputs either haplotypes or founder sequences, one sequence per file. The files will be written to the current working directory.

### Each output sequence requires a file descriptor per one pass of the variant graph

When generating sample sequences, one sequence is generated per each sample and chromosome copy. Similarly, when generating founder sequences, one seuqence is generated per each founder. To this end, `variant_graph_to_sequences` attempts to open as many files before parsing the variant file. If the number of samples or founders is high, the number of files may exceed shell or operating system limits. In order to process the variant graph as few times as possible, the maximum number of file descriptors available to processes started by the shell may be changed with `ulimit -n`. For example, `ulimit -n 8192` sets the maximum number to 8192. Another option is to reduce the number of samples to be handled in one pass with vcf2multialign’s `--chunk-size` command line option.

### Inspecting the variant graph

Two additional tools are provided for variant graph inspection. `variant-graph-to-gv/variant_graph_to_gv` produces [Graphviz](https://www.graphviz.org) output from a variant graph. For anything except very small graphs, `inspect-variant-graph/inspect_variant_graph` may be more suitable as it outputs only one node at a time.

## Creating unaligned sample sequences

Instead of generating a reference-guided multiple alignment, unaligned sequences may be generated with `vcf-to-unaligned/vcf_to_unaligned` . The sample and chromosome copy numbers need to be specified as command line parameters.

## Creating a variant call file from a pair of aligned sequences and a variant file

`combine-msa-vcf/combine_msa_vcf` may be used to generate a variant call file from two aligned sequences. Optionally, variants relative to the second sequence may be included in the resulting file; they will be merged with the variants generated from the alignment and rewritten to be relative to the first sequence.
