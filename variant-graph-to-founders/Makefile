include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				sequence_generator_base.o

all: variant_graph_to_founders

clean:
	$(RM) $(OBJECTS) variant_graph_to_founders cmdline.c cmdline.h version.h

variant_graph_to_founders: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
