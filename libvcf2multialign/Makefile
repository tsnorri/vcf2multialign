include ../local.mk
include ../common.mk

OBJECTS =	variant_graph.o \
			transpose_matrix.o


all: libvcf2multialign.a

clean:
	$(RM) $(OBJECTS) libvcf2multialign.a
	
libvcf2multialign.a: $(OBJECTS)
	$(AR) rcs $@ $^
