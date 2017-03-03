# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

#CC			?= gcc
#CXX			?= g++
CC			= clang-3.7
CXX			= clang++-3.7
GENGETOPT	?= gengetopt

CPPFLAGS	= -I../include -I../lib/vcflib/include -I../lib/vcflib/tabixpp/htslib -I/sw/opt/boost-1_58/include
LDFLAGS		= -L/sw/opt/boost-1_58/lib -lboost_iostreams -L../lib/vcflib -lvcflib -L../lib/vcflib/tabixpp -ltabix++ -L../lib/vcflib/tabixpp/htslib -lhts -lz
CXXFLAGS	= -O0 -g -std=c++14
CFLAGS		= -O0 -g -std=c99

OBJECTS		=	cmdline.o \
				main.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS)

vcf2multialign: $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS)


%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"


main.c : cmdline.c
