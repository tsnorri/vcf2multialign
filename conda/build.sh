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
#cxx_includes="-isystem '${PREFIX}/include/c++/v1'"
clang_version=`${PREFIX}/bin/clang -dumpversion`
includes="-I '${SRC_DIR}/conda'"
cxx_includes="-isystem '${PREFIX}/x86_64-conda-linux-gnu/include/c++/9.3.0' -isystem '${PREFIX}/x86_64-conda-linux-gnu/include/c++/9.3.0/x86_64-conda-linux-gnu'"
sys_includes="-isystem '${PREFIX}/lib/clang/${clang_version}/include' -isystem '${PREFIX}/x86_64-conda-linux-gnu/sysroot/usr/include'"
cppflags="-nostdinc -U__STDC_HOSTED__ ${includes}"
system_cppflags="${cxx_includes} ${sys_includes}"
cflags="-fblocks"
cxxflags="-fblocks -nostdinc++"

echo "CC = ${clang}"																> local.mk
echo "CXX = ${clangxx}"																>> local.mk
echo "CPPFLAGS = ${cppflags}"														>> local.mk
echo "CFLAGS = ${cflags}"															>> local.mk
echo "CXXFLAGS = ${cxxflags}"														>> local.mk
echo "LDFLAGS = -L${PREFIX}/lib -ldl"												>> local.mk
echo "SYSTEM_CPPFLAGS = ${system_cppflags}"											>> local.mk
echo "LIBDISPATCH_CFLAGS = ${cppflags} ${sys_includes}"								>> local.mk
echo "LIBDISPATCH_CXXFLAGS = ${cppflags} -nostdinc ${cxx_includes} ${sys_includes}"	>> local.mk
echo "BOOST_ROOT = ${PREFIX}"														>> local.mk
echo "BOOST_LIBS = -lboost_iostreams"												>> local.mk
echo "LIBBSD_LIB = ${PREFIX}/x86_64-conda-linux-gnu/sysroot/usr/lib64/libbsd.a"		>> local.mk

make -j ${CPU_COUNT}

dst_bin="${PREFIX}/bin"
mkdir -p "${dst_bin}"
cp combine-msa-vcf/combine_msa_vcf 							"${dst_bin}"
cp create-variant-graph/create_variant_graph				"${dst_bin}"
cp preprocess-vcf/preprocess_vcf							"${dst_bin}"
cp variant-graph-to-sequences/variant_graph_to_sequences	"${dst_bin}"
cp vcf-to-unaligned/vcf_to_unaligned						"${dst_bin}"
