include ../local.mk
include ../common.mk

CATCH2_PREFIX			= ../lib/libbio/lib/Catch2
RAPIDCHECK_PREFIX		= ../lib/libbio/lib/rapidcheck

CPPFLAGS +=	-I../lib/libbio/lib/Catch2/src \
            -I../lib/libbio/lib/Catch2/build/generated-includes \
            -I../lib/libbio/lib/rapidcheck/include \
            -I../lib/libbio/lib/rapidcheck/extras/catch/include

OBJECTS	=	founder_sequences.o \
			transpose_matrix.o \
			variant_graph.o \
			main.o

ifeq ($(origin NO_COVERAGE_CHECK),undefined)
GCOV_LIBRARY = -lgcov
LIBVCF2MULTIALIGN_PATH ?= ../libvcf2multialign/libvcf2multialign.coverage.a
else
GCOV_LIBRARY =
LIBVCF2MULTIALIGN_PATH ?= ../libvcf2multialign/libvcf2multialign.a
endif

DEPENDENCY_LIBRARIES = $(LIBVCF2MULTIALIGN_PATH) ../lib/libbio/src/libbio.a $(RAPIDCHECK_PREFIX)/build/librapidcheck.a $(CATCH2_PREFIX)/build/src/libCatch2Main.a $(CATCH2_PREFIX)/build/src/libCatch2.a
TEST_LDFLAGS = $(GCOV_LIBRARY) $(BOOST_LIBS) $(LDFLAGS)


all: tests


.PHONY: run-tests clean clean-tests


run-tests: tests
	./tests


coverage: coverage/index.html


tests: $(OBJECTS) $(DEPENDENCY_LIBRARIES)
	$(CXX) -o $@ $^ $(TEST_LDFLAGS)


clean: clean-tests
	$(RM) $(OBJECTS) tests
	$(RM) -r coverage


coverage/index.html: run-tests
	$(MKDIR) -p coverage
	$(GCOVR) --gcov-executable $(GCOV) --root .. --filter ../include --filter ../libvcf2multialign --print-summary --html --html-details -o coverage/index.html ../libvcf2multialign
