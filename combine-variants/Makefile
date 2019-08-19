include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				path_sorted_variant.o \
				preprocess_variants.o \
				sample_sorter.o \
				variant_preprocessor.o

all: combine_variants

clean:
	$(RM) $(OBJECTS) combine_variants cmdline.c cmdline.h version.h vcf_reader.pdf

combine_variants: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
