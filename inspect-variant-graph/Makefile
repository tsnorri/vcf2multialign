include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: inspect_variant_graph

clean:
	$(RM) $(OBJECTS) inspect_variant_graph cmdline.c cmdline.h version.h

inspect_variant_graph: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
