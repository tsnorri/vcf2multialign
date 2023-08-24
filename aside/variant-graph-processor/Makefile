include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: variant_graph_processor

clean:
	$(RM) $(OBJECTS) variant_graph_processor cmdline.c cmdline.h version.h config.h

variant_graph_processor: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
