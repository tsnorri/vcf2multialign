include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				path_sorted_variant.o \
				preprocess_variants.o \
				sample_sorter.o \
				variant_preprocessor.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h vcf_reader.pdf

vcf2multialign: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../libvcf2multialign/libvcf2multialign.a

main.c : cmdline.c
cmdline.c : config.h

include ../config.mk
