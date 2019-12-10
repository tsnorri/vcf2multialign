include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: index_vcf

clean:
	$(RM) $(OBJECTS) index_vcf cmdline.c cmdline.h version.h

index_vcf: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
