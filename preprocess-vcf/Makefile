include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				preprocess_variants.o

all: preprocess_vcf

clean:
	$(RM) $(OBJECTS) preprocess_vcf cmdline.c cmdline.h version.h

preprocess_vcf: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
