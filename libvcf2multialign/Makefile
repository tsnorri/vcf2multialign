include ../local.mk
include ../common.mk

OBJECTS = 	check_overlapping_non_nested_variants.o \
			generate_context_base.o \
			read_single_fasta_seq.o \
			variant_handler_base.o


all: libvcf2multialign.a

clean:
	$(RM) $(OBJECTS) libvcf2multialign.a
	
libvcf2multialign.a: $(OBJECTS)
	$(AR) rcs $@ $^
