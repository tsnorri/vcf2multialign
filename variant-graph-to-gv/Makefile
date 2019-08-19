include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: variant_graph_to_gv

clean:
	$(RM) $(OBJECTS) variant_graph_to_gv cmdline.c cmdline.h version.h

variant_graph_to_gv: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
