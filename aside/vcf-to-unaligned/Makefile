include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: vcf_to_unaligned

clean:
	$(RM) $(OBJECTS) vcf_to_unaligned cmdline.c cmdline.h version.h

vcf_to_unaligned: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
