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
