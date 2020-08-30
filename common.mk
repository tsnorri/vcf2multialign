# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

# Default values.
WARNING_FLAGS	?= -Wall -Werror -Wno-deprecated-declarations -Wno-unused
OPT_FLAGS		?= -O2 -g

CMAKE			?= cmake
CP				?= cp
DOT				?= dot
GENGETOPT		?= gengetopt
MKDIR			?= mkdir
NINJA			?= ninja
RAGEL			?= ragel
TAR				?= tar
WGET			?= wget

CFLAGS			?=
CXXFLAGS		?=
CPPFLAGS		?=
LDFLAGS			?=
SYSTEM_CFLAGS	?=
SYSTEM_CXXFLAGS	?=
SYSTEM_CPPFLAGS	?=
SYSTEM_LDFLAGS	?=

# Target type description, used currently in the .tar.gz name.
TARGET_TYPE		?=

BOOST_ROOT		?= /usr
BOOST_INCLUDE	?= -I$(BOOST_ROOT)/include

CFLAGS			+= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CFLAGS)
CXXFLAGS		+= -std=c++2a $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CXXFLAGS)
CPPFLAGS		+= $(SYSTEM_CPPFLAGS) -DHAVE_CONFIG_H -I../include -I../lib/cereal/include -I../lib/libbio/include -I../lib/libbio/lib/GSL/include -I../lib/libbio/lib/range-v3/include $(BOOST_INCLUDE)
LDFLAGS			+= $(SYSTEM_LDFLAGS) ../lib/libbio/src/libbio.a $(BOOST_LIBS)

ifeq ($(shell uname -s),Linux)
	CPPFLAGS	+= -I../lib/swift-corelibs-libdispatch
	LDFLAGS		+= ../lib/swift-corelibs-libdispatch/build/src/libdispatch.a ../lib/swift-corelibs-libdispatch/build/src/BlocksRuntime/libBlocksRuntime.a /usr/lib/x86_64-linux-gnu/libbsd.a -lpthread -lz
endif


# FIXME: the first two likely only work with Clang; I think GCC uses something else than -coverage.
%.cov.o: %.c
	$(CC) -c -coverage $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.cov.o: %.cc
	$(CXX) -c -coverage $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"

%.cc: %.rl
	$(RAGEL) -L -C -G2 -o $@ $<

%.dot: %.rl
	$(RAGEL) -V -p -o $@ $<

%.pdf: %.dot
	$(DOT) -Tpdf $< > $@
