include ../local.mk
include ../common.mk

CATCH2_PREFIX			= ../lib/libbio/lib/Catch2
CATCH2_HEADERS			= $(shell find $(CATCH2_PREFIX)/include)

CPPFLAGS += -I$(CATCH2_PREFIX)/single_include
CXXFLAGS += -coverage
LDFLAGS += -coverage

OBJECTS	=	founder_sequences.o \
			transpose_matrix.o \
			variant_graph.o \
			main.o

GCDA = $(OBJECTS:.o=.gcda)
GCNO = $(OBJECTS:.o=.gcno)


LIBRARIES = ../libvcf2multialign/libvcf2multialign.coverage.a ../lib/libbio/src/libbio.a ../lib/libbio/lib/rapidcheck/build/librapidcheck.a


all: tests


.PHONY: run-tests clean clean-tests coverage


run-tests: tests
	./tests


coverage: run-tests
	$(MKDIR) -p coverage
	$(GCOVR) --gcov-executable $(GCOV) --root .. --filter ../include --filter ../libvcf2multialign --print-summary --html --html-details -o coverage/index.html . ../libvcf2multialign


tests: $(OBJECTS) $(LIBRARIES)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBRARIES) $(BOOST_LIBS) -lgcov


clean: clean-tests
	$(RM) $(OBJECTS) $(GCDA) tests


clean-tests:
	$(RM) $(GCNO) test-coverage.covdata
	$(RM) -r test-coverage


$(CATCH2_PREFIX)/single_include/catch2/catch.hpp: $(CATCH2_HEADERS)
	cd $(CATCH2_PREFIX) && $(PYTHON) scripts/generateSingleHeader.py
