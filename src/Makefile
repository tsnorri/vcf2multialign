include ../local.mk
include ../common.mk

OBJECTS		=	check_overlapping_non_nested_variants.o \
				cmdline.o \
				error_logger.o \
				generate_context_base.o \
				generate_graph_context.o \
				generate_haplotypes.o \
				generate_sequences_context.o \
				graph_writer.o \
				main.o \
				read_single_fasta_seq.o \
				variant_buffer.o \
				variant_handler.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h vcf_reader.pdf

vcf2multialign: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)

main.c : cmdline.c
cmdline.c : config.h

include ../config.mk
