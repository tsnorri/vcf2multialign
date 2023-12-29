#!/bin/bash

set -e
set -x

if [ -f local.mk ]
then
	echo "ERROR: local.mk already exists."
	exit 1
fi

# FIXME: hack to make conda-build use Clang on Linux. This should be made platform-independent.
clang="${PREFIX}/bin/clang"
clangxx="${PREFIX}/bin/clang++"
clang_version=`${PREFIX}/bin/clang -dumpversion`
includes="--sysroot '${PREFIX}/x86_64-conda-linux-gnu/sysroot' -isystem =/../../include/c++/13.2.0 -isystem =/../include -isystem =/../../include" # -I'${SRC_DIR}/conda'"
cppflags="${includes}"
cflags="-fblocks -Wno-unknown-warning-option"
cxxflags="-fblocks -Wno-unknown-warning-option"

echo "CC = ${clang}"																> local.mk
echo "CXX = ${clangxx}"																>> local.mk
echo "CPPFLAGS = ${cppflags}"														>> local.mk
echo "CFLAGS = ${cflags}"															>> local.mk
echo "CXXFLAGS = ${cxxflags}"														>> local.mk
echo "LDFLAGS = -L${PREFIX}/lib -ldl"												>> local.mk
echo "SYSTEM_CPPFLAGS = ${cppflags}"												>> local.mk
#echo "LIBDISPATCH_CFLAGS = ${cppflags}"												>> local.mk
#echo "LIBDISPATCH_CXXFLAGS = ${cppflags}"											>> local.mk
echo "BOOST_ROOT = ${PREFIX}"														>> local.mk
echo "BOOST_INCLUDE ="																>> local.mk
echo "BOOST_LIBS = -lboost_iostreams"												>> local.mk
echo "LIBBSD_LIB = ${PREFIX}/x86_64-conda-linux-gnu/sysroot/usr/lib64/libbsd.a"		>> local.mk

cat local.mk

make #-j ${CPU_COUNT}

dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp combine-msa-vcf/combine_msa_vcf 							"${dst_bin}"
cp create-variant-graph/create_variant_graph				"${dst_bin}"
cp preprocess-vcf/preprocess_vcf							"${dst_bin}"
cp variant-graph-to-sequences/variant_graph_to_sequences	"${dst_bin}"
cp vcf-to-unaligned/vcf_to_unaligned						"${dst_bin}"
