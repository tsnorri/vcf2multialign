include ../local.mk
include ../common.mk

OBJECTS =	can_handle_variant_alts.o \
			check_ploidy.o \
			log_assertion_failure.o \
			path_sorted_variant.o \
			progress_indicator_manager.o \
			read_single_fasta_seq.o \
			sample_sorter.o \
			variant_format.o \
			variant_graph.o \
			variant_graph_generator.o \
			variant_partitioner.o \
			variant_processor_delegate.o \
			vcf_processor.o


all: libvcf2multialign.a

clean:
	$(RM) $(OBJECTS) libvcf2multialign.a
	
libvcf2multialign.a: $(OBJECTS)
	$(AR) rcs $@ $^
