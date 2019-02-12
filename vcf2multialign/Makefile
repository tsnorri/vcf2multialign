include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				error_logger.o \
				generate_graph_context.o \
				generate_haplotypes.o \
				generate_sequences_context.o \
				graph_writer.o \
				main.o \
				variant_buffer.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h vcf_reader.pdf

vcf2multialign: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../libvcf2multialign/libvcf2multialign.a

main.c : cmdline.c
cmdline.c : config.h

include ../config.mk
