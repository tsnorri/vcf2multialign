include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h config.h

vcf2multialign: $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
