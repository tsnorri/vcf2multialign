include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				find_optimal_cut_positions.o \
				main.o

all: find_optimal_cut_positions

clean:
	$(RM) $(OBJECTS) find_optimal_cut_positions cmdline.c cmdline.h version.h

find_optimal_cut_positions: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
