include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				create_variant_graph.o \
				main.o

all: create_variant_graph

clean:
	$(RM) $(OBJECTS) create_variant_graph cmdline.c cmdline.h version.h config.h

create_variant_graph: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
