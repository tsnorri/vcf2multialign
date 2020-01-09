include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				combine_msa.o \
				main.o \
				msa_combiner.o \
				msa_combiner_parse_msa.o \
				vcf_record_generator.o

all: combine_msa_vcf

clean:
	$(RM) $(OBJECTS) combine_msa_vcf cmdline.c cmdline.h version.h

combine_msa_vcf: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
