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
make -j ${CPU_COUNT} all lib/libbio/vcfcat/vcfcat

echo "Copying build products"
dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp vcf2multialign/vcf2multialign	"${dst_bin}"
cp lib/libbio/vcfcat/vcfcat			"${dst_bin}"

echo "Copying documentation"
dst_doc="${PREFIX}/share/doc/vcf2multialign"
mkdir -p "${dst_doc}"
cp README.md "${dst_doc}"
cp LICENSE "${dst_doc}"
cp lib/cereal/LICENSE "${dst_doc}/cereal-license.txt"
cp lib/libbio/lib/GSL/LICENSE "${dst_doc}/GSL-license.txt"
cp lib/libbio/lib/range-v3/LICENSE.txt "${dst_doc}/range-v3-license.txt"
