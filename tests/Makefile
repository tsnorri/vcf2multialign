include ../local.mk
include ../common.mk

CATCH2_PREFIX			= ../lib/libbio/lib/Catch2
CATCH2_HEADERS			= $(shell find $(CATCH2_PREFIX)/include)

CPPFLAGS += -I$(CATCH2_PREFIX)/single_include

OBJECTS	=	founder_sequences.o \
			transpose_matrix.o \
			variant_graph.o \
			main.o


LIBRARIES = ../libvcf2multialign/libvcf2multialign.coverage.a ../lib/libbio/src/libbio.a ../lib/libbio/lib/rapidcheck/build/librapidcheck.a


all: tests


.PHONY: run-tests clean clean-tests


run-tests: tests
	./tests


coverage: coverage/index.html


tests: $(OBJECTS) $(LIBRARIES)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBRARIES) $(BOOST_LIBS) -lgcov


clean: clean-tests
	$(RM) $(OBJECTS) tests
	$(RM) -r coverage


$(CATCH2_PREFIX)/single_include/catch2/catch.hpp: $(CATCH2_HEADERS)
	cd $(CATCH2_PREFIX) && $(PYTHON) scripts/generateSingleHeader.py


coverage/index.html: run-tests
	$(MKDIR) -p coverage
	$(GCOVR) --gcov-executable $(GCOV) --root .. --filter ../include --filter ../libvcf2multialign --print-summary --html --html-details -o coverage/index.html ../libvcf2multialign
