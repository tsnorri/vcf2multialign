include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: variant_graph_statistics

clean:
	$(RM) $(OBJECTS) variant_graph_statistics cmdline.c cmdline.h version.h config.h

variant_graph_statistics: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
