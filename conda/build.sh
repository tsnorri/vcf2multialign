#!/bin/bash

set -euxo pipefail

if [ -f local.mk ]
then
	echo "ERROR: local.mk already exists."
	exit 1
fi

echo "Generating local.mk"
m4 -D CONDA_PREFIX="${PREFIX}" conda/local.mk.m4 > local.mk

echo "Running make"
make -j ${CPU_COUNT}

echo "Copying build products"
dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp combine-msa-vcf/combine_msa_vcf 							"${dst_bin}"
cp create-variant-graph/create_variant_graph				"${dst_bin}"
cp preprocess-vcf/preprocess_vcf							"${dst_bin}"
cp variant-graph-to-sequences/variant_graph_to_sequences	"${dst_bin}"
cp vcf-to-unaligned/vcf_to_unaligned						"${dst_bin}"
