include ../local.mk
include ../common.mk

.PRECIOUS: vcf_reader.cc

OBJECTS		=	all_haplotypes_task.o \
				check_overlapping_non_nested_variants.o \
				cmdline.o \
				error_logger.o \
				file_handling.o \
				find_subgraph_starting_points.o \
				generate_haplotypes.o \
				main.o \
				preparation_task.o \
				progress_bar.o \
				read_single_fasta_seq.o \
				sequence_writer.o \
				status_logger.o \
				types.o \
				util.o \
				variant_buffer.o \
				variant_handler.o \
				variant_stats.o \
				variant.o \
				vcf_input.o \
				vcf_reader.o

all: vcf2multialign

clean:
	$(RM) $(OBJECTS) vcf2multialign cmdline.c cmdline.h version.h vcf_reader.pdf

vcf2multialign: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS)

test: test.o
	$(CXX) -o test test.o $(LDFLAGS)

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
