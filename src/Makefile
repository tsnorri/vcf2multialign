include ../local.mk
include ../common.mk

.PRECIOUS: vcf_reader.cc

OBJECTS		=	check_overlapping_non_nested_variants.o \
				cmdline.o \
				error_logger.o \
				generate_haplotypes.o \
				main.o \
				read_single_fasta_seq.o \
				variant_buffer.o \
				variant_handler.o \
				variant.o \
				vcf_reader.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS)

vcf2multialign: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)

main.c : cmdline.c
