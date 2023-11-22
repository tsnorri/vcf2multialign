include ../local.mk
include ../common.mk

OBJECTS =	find_cut_positions.o \
			founder_sequence_greedy_output.o \
			haplotype_output.o \
			output.o \
			sequence_writer.o \
			transpose_matrix.o \
			variant_graph.o

OBJECTS_COVERAGE	= $(OBJECTS:.o=.cov.o)
GCDA				= $(OBJECTS_COVERAGE:.o=.gcda)
GCNO				= $(OBJECTS_COVERAGE:.o=.gcno)


all: libvcf2multialign.a

clean:
	$(RM) $(OBJECTS) $(OBJECTS_COVERAGE) $(GCDA) $(GCNO) libvcf2multialign.a libvcf2multialign.coverage.a
	
libvcf2multialign.a: $(OBJECTS)
	$(AR) rcs $@ $^

libvcf2multialign.coverage.a: $(OBJECTS_COVERAGE)
	$(AR) rcs $@ $^
