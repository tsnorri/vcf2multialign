# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

# Default values.
# Boost uses some deprecated builtins (as of Clang 14).
WARNING_FLAGS	?= -Wall -Werror -Wno-deprecated-declarations -Wno-unused
OPT_FLAGS		?= -O2 -g

CMAKE			?= cmake
CP				?= cp
DOT				?= dot
GENGETOPT		?= gengetopt
MKDIR			?= mkdir
NINJA			?= ninja
PATCH			?= patch
RAGEL			?= ragel
TAR				?= tar
WGET			?= wget
GCOV			?= gcov
GCOVR			?= gcovr

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
BOOST_LIBS		?= -lboost_iostreams

CFLAGS			+= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CFLAGS)
CXXFLAGS		+= -std=c++2b $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CXXFLAGS)
CPPFLAGS		+= -DHAVE_CONFIG_H -I../include -I../lib/cereal/include -I../lib/libbio/include -I../lib/libbio/lib/GSL/include -I../lib/libbio/lib/range-v3/include $(BOOST_INCLUDE) $(SYSTEM_CPPFLAGS)
LDFLAGS			:= $(LDFLAGS) $(SYSTEM_LDFLAGS)

%.cov.o: %.c
	$(CC) -c --coverage $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.cov.o: %.cc
	$(CXX) -c --coverage $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

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
