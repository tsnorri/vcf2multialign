# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

WARNING_FLAGS	= -Wall -Werror -Wno-deprecated-declarations -Wno-unused
OPT_FLAGS		= -O2 -g

CFLAGS			= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS)
CXXFLAGS		= -std=c++17 $(OPT_FLAGS) $(WARNING_FLAGS)
CPPFLAGS		= -DHAVE_CONFIG_H -I../include -I../lib/libbio/include -I../lib/libdispatch -I../lib/libpwq/include $(BOOST_INCLUDE)
LDFLAGS			= ../lib/libbio/src/libbio.a $(LIBDISPATCH_LIBS) $(BOOST_LIBS)

CMAKE			?= cmake
CP				?= cp
DOT				?= dot
GENGETOPT		?= gengetopt
MKDIR			?= mkdir
RAGEL			?= ragel
RM				?= rm


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
