# vcf2multialign

Create a set of reference-guided multiple-aligned haplotypes or founder sequences from a variant call file and a reference sequence. Please see [releases](https://github.com/tsnorri/vcf2multialign/releases) and [Anaconda.org](https://anaconda.org/tsnorri/vcf2multialign) for pre-built binaries.

## Academic Use

If you use the software in an academic setting we kindly ask you to cite [Founder reconstruction enables scalable and seamless pangenomic analysis](https://doi.org/10.1093/bioinformatics/btab516).

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

## Building

To clone the repository with submodules, please use `git clone --recursive https://github.com/tsnorri/vcf2multialign.git`.

### With [conda-build](https://docs.conda.io/projects/conda-build/en/stable/index.html)

A conda package can be built with conda-build as follows. The build script has been tested with conda-build 3.27.0. [glibc](https://www.gnu.org/software/libc/) 2.28 or newer is required.

1. `cd conda`
2. `./conda-build.sh`

### By hand

The following software and libraries are required to build vcf2multialign. The tested versions are also listed.

- Reasonably new compilers for C and C++. We use [GCC 12.3](https://gcc.gnu.org) in C++2b mode for our builds.
- [GNU gengetopt 2.23](https://www.gnu.org/software/gengetopt/gengetopt.html) (tested with version 2.22.6)
- [Ragel State Machine Compiler 6.10](http://www.colm.net/open-source/ragel/)
- [Boost 1.82.0](http://www.boost.org)
- [libbsd](https://libbsd.freedesktop.org/) on Linux.

After installing the prerequisites, please do the following:

1. Create a file called `local.mk` in the root of the cloned repository to specify build variables. One of the files [linux-static.local.mk](linux-static.local.mk) and [conda/local.mk.m4](conda/local.mk.m4) may be used as a starting point.
2. Run Make with e.g. `make -j16 dist` to create a gzipped tar archive of the executables.

## Usage

Outputting a reference-guided multiple sequence alignment of predicted haplotype sequences to `haplotypes.a2m` (with only uppercase characters). Sequence `1` from `hs37d5.fa` is used as the reference and the variants of the `chr1` chromosome from `variants.vcf` as the variants:

```
vcf2multialign --haplotypes --input-reference=hs37d5.fa --reference-sequence=1 --input-variants=variants.vcf --output-sequences-a2m=founders.a2m --chromosome=chr1
```

Outputting a reference-guided multiple sequence alignment of 25 founder sequences to `founders.a2m` with a minimum aligned distance of 50 between the graph components using the same inputs:

```
vcf2multialign --founder-sequences=25 --minimum-distance=50 --input-reference=hs37d5.fa --reference-sequence=1 --input-variants=variants.vcf --output-sequences-a2m=founders.a2m --chromosome=chr1
```

Please refer to `vcf2multialign --help` for a complete list of options.
