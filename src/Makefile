include ../local.mk
include ../common.mk

OBJECTS		=	check_overlapping_non_nested_variants.o \
				cmdline.o \
				error_logger.o \
				generate_haplotypes.o \
				main.o \
				read_single_fasta_seq.o \
				variant_buffer.o \
				variant_handler.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h vcf_reader.pdf

vcf2multialign: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)

main.c : cmdline.c
cmdline.c : config.h

# Use the name config.h since cmdline.c includes it automatically.
config.h: SHELL := /bin/bash
config.h: Makefile ../.git
	IFS=$$'\n'; \
	TAGS=`git tag -l --points-at HEAD`; \
	HASH=`git rev-parse --short --verify HEAD`; \
	VERSION="DEV"; \
	if [ -n "$${TAGS}" ]; \
	then \
		for name in $${TAGS}; \
		do \
			if [[ "$${name}" == v* ]]; \
			then \
				VERSION="$${name:1}"; \
				break; \
			fi; \
		done; \
	fi; \
	printf "#define CMDLINE_PARSER_VERSION \"%s %s\"\n" "$${VERSION}" "$${HASH}" > config.h
