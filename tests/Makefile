include ../local.mk
include ../common.mk

CATCH2_PREFIX			= ../lib/libbio/lib/Catch2
CATCH2_HEADERS			= $(shell find $(CATCH2_PREFIX)/include)

CPPFLAGS += -I$(CATCH2_PREFIX)/single_include
CXXFLAGS += -coverage
LDFLAGS += -coverage

OBJECTS	=	combine_msa_vcf.o \
			main.o

GCDA = $(OBJECTS:.o=.gcda)
GCNO = $(OBJECTS:.o=.gcno)


LIBRARIES = ../lib/libbio/src/libbio.a ../libvcf2multialign/libvcf2multialign.a ../combine-msa-vcf/libcombinemsa.coverage.a


all: run-tests


.PHONY: run-tests clean clean-tests


run-tests: tests
	./tests
	lcov --capture --directory . --directory ../combine-msa-vcf --output-file test-coverage.covdata --include '*/combine-msa-vcf/*'
	genhtml --legend --title test-coverage -o test-coverage test-coverage.covdata


tests: $(OBJECTS) $(LIBRARIES)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBRARIES) $(BOOST_LIBS)


clean: clean-tests
	$(RM) $(OBJECTS) $(GCDA) tests


clean-tests:
	$(RM) $(GCNO) test-coverage.covdata
	$(RM) -r test-coverage


$(CATCH2_PREFIX)/single_include/catch2/catch.hpp: $(CATCH2_HEADERS)
	cd $(CATCH2_PREFIX) && $(PYTHON) scripts/generateSingleHeader.py
