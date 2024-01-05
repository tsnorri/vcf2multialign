#!/bin/bash

set -e
set -x

if [ -f local.mk ]
then
	echo "ERROR: local.mk already exists."
	exit 1
fi

echo "Generating local.mk"
m4 -D CONDA_PREFIX="${PREFIX}" conda/local.mk.m4 > local.mk

echo "Contents of local.mk:"
cat local.mk

pushd "${PREFIX}/x86_64-conda-linux-gnu/sysroot/usr/lib64/"
ln -s ../../../../lib/gcc/x86_64-conda-linux-gnu/13.2.0/crtbeginS.o ./
ln -s ../../../../lib/gcc/x86_64-conda-linux-gnu/13.2.0/crtendS.o ./
popd

make #-j ${CPU_COUNT}

dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp combine-msa-vcf/combine_msa_vcf 							"${dst_bin}"
cp create-variant-graph/create_variant_graph				"${dst_bin}"
cp preprocess-vcf/preprocess_vcf							"${dst_bin}"
cp variant-graph-to-sequences/variant_graph_to_sequences	"${dst_bin}"
cp vcf-to-unaligned/vcf_to_unaligned						"${dst_bin}"
