include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				path_sorted_variant.o \
				preprocess_variants.o \
				sample_sorter.o \
				variant_preprocessor.o

all: preprocess_vcf

clean:
	$(RM) $(OBJECTS) preprocess_vcf cmdline.c cmdline.h version.h vcf_reader.pdf

preprocess_vcf: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
