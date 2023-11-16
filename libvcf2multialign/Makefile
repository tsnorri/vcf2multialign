include ../local.mk
include ../common.mk

OBJECTS =	find_cut_positions.o \
			founder_sequence_greedy_output.o \
			haplotype_output.o \
			output.o \
			sequence_writer.o \
			transpose_matrix.o \
			variant_graph.o


all: libvcf2multialign.a

clean:
	$(RM) $(OBJECTS) libvcf2multialign.a
	
libvcf2multialign.a: $(OBJECTS)
	$(AR) rcs $@ $^
