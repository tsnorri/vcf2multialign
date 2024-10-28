include ../local.mk
include ../common.mk

#CPPFLAGS += -DLIBBIO_LOG_ALLOCATED_MEMORY

OBJECTS		=	cmdline.o \
				log_memory_usage.o \
				main.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h config.h

vcf2multialign: $(OBJECTS) ../libvcf2multialign/libvcf2multialign.a ../lib/libbio/src/libbio.a
	$(CXX) -o $@ $^ $(BOOST_LIBS) $(LDFLAGS)

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
