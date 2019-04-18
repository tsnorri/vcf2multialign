include ../local.mk
include ../common.mk

OBJECTS =	can_handle_variant_alts.o \
			check_ploidy.o \
			read_single_fasta_seq.o \
			variant_format.o


all: libvcf2multialign.a

clean:
	$(RM) $(OBJECTS) libvcf2multialign.a
	
libvcf2multialign.a: $(OBJECTS)
	$(AR) rcs $@ $^
