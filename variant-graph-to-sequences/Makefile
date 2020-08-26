include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				founder_sequence_greedy_generator.o \
				main.o \
				sequence_generator_base.o

all: variant_graph_to_sequences

clean:
	$(RM) $(OBJECTS) variant_graph_to_sequences cmdline.c cmdline.h version.h

variant_graph_to_sequences: $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
