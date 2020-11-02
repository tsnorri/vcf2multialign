include ../local.mk
include ../common.mk

OBJECTS				=	cmdline.o \
						main.o
						
LIBRARY_OBJECTS		=	combine_msa.o \
						mnv_combiner.o \
						msa_combiner.o \
						msa_data_source.o \
						msa_data_source_parse_msa.o \
						overlap_counter.o \
						utility.o \
						variant_filter.o \
						variant_writer.o \
						vcf_record_generator.o

LIBRARY_OBJECTS_COV	=	$(LIBRARY_OBJECTS:.o=.cov.o)
GCDA				=	$(LIBRARY_OBJECTS_COV:.o=.gcda)
GCNO				=	$(LIBRARY_OBJECTS_COV:.o=.gcno)


all: combine_msa_vcf
	
.PHONY: clean

.PRECIOUS: msa_data_source_parse_msa.cc

clean:
	$(RM)	$(OBJECTS) \
			$(LIBRARY_OBJECTS) \
			$(LIBRARY_OBJECTS_COV) \
			$(GCDA) \
			$(GCNO) \
			combine_msa_vcf \
			libcombinemsa.a \
			libcombinemsa.coverage.a \
			msa_data_source_parse_msa.cc \
			cmdline.c \
			cmdline.h \
			version.h

combine_msa_vcf: $(OBJECTS) libcombinemsa.a
	$(CXX) -o $@ $(OBJECTS) libcombinemsa.a ../libvcf2multialign/libvcf2multialign.a $(LDFLAGS)

libcombinemsa.a: $(LIBRARY_OBJECTS)
	$(AR) rcs $@ $^

libcombinemsa.coverage.a: $(LIBRARY_OBJECTS_COV)
	$(AR) rcs $@ $^


main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
