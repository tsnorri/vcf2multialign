include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				vcf_record_generator.o

all: vcfcontains

clean:
	$(RM) $(OBJECTS) vcfcontains cmdline.c cmdline.h version.h

vcfcontains: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
