# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

# Default values.
WARNING_FLAGS	?= -Wall -Werror -Wno-deprecated-declarations -Wno-unused
OPT_FLAGS		?= -O2 -g

CFLAGS			?=
CXXFLAGS		?=
CPPFLAGS		?=
LDFLAGS			?=
SYSTEM_CFLAGS	?=
SYSTEM_CXXFLAGS	?=
SYSTEM_CPPFLAGS	?=
SYSTEM_LDFLAGS	?=

AR				?= ar
CC				?= cc
CMAKE			?= cmake
CP				?= cp
CXX				?= c++
DOT				?= dot
GENGETOPT		?= gengetopt
MKDIR			?= mkdir
RAGEL			?= ragel
RM				?= rm

BOOST_INCLUDE	?= /usr/include

CFLAGS			+= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CFLAGS)
CXXFLAGS		+= -std=c++17 $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CXXFLAGS)
CPPFLAGS		+= -DHAVE_CONFIG_H -I../include -I../lib/libbio/include -I../lib/libbio/lib/GSL/include -I../lib/libbio/lib/range-v3/include -I../lib/msa2dag/include -I../lib/libdispatch -I../lib/libpwq/include $(BOOST_INCLUDE)
LDFLAGS			+= ../lib/libbio/src/libbio.a ../lib/msa2dag/lib/libMsa2Dag.a $(LIBDISPATCH_LIBS) $(BOOST_LIBS)


%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"

%.cc: %.rl
	$(RAGEL) -L -C -G2 -o $@ $<

%.dot: %.rl
	$(RAGEL) -V -p -o $@ $<

%.pdf: %.dot
	$(DOT) -Tpdf $< > $@
