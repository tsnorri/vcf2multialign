include local.mk
include common.mk

DEPENDENCIES = lib/msa2dag/lib/libMsa2Dag.a lib/libbio/src/libbio.a
ifeq ($(shell uname -s),Linux)
	DEPENDENCIES += lib/swift-corelibs-libdispatch/build/src/libdispatch.a
endif

# “$() $()” is a literal space.
OS_NAME = $(shell tools/os_name.sh)
VERSION = $(subst $() $(),-,$(shell tools/git_version.sh))
DIST_TARGET_DIR = vcf2multialign-$(VERSION)
DIST_NAME_SUFFIX = $(if $(TARGET_TYPE),-$(TARGET_TYPE),)
DIST_TAR_GZ = vcf2multialign-$(VERSION)-$(OS_NAME)$(DIST_NAME_SUFFIX).tar.gz


.PHONY: all clean-all clean clean-dependencies dependencies

all: libvcf2multialign/libvcf2multialign.a preprocess-vcf/preprocess_vcf vcf2multialign/vcf2multialign

clean-all: clean clean-dependencies clean-dist

clean:
	$(MAKE) -C libvcf2multialign clean
	$(MAKE) -C preprocess-vcf clean
	$(MAKE) -C vcf2multialign clean

clean-dependencies: lib/libbio/local.mk
	$(MAKE) -C lib/libbio clean-all
	$(MAKE) -C lib/msa2dag clean-all
	$(RM) -r lib/swift-corelibs-libdispatch/build

clean-dist:
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

dist: $(DIST_TAR_GZ)

preprocess-vcf/preprocess_vcf: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C preprocess-vcf

vcf2multialign/vcf2multialign: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C vcf2multialign

libvcf2multialign/libvcf2multialign.a: $(DEPENDENCIES)
	$(MAKE) -C libvcf2multialign

$(DIST_TAR_GZ): preprocess-vcf/preprocess_vcf vcf2multialign/vcf2multialign
	$(MKDIR) -p $(DIST_TARGET_DIR)
	$(CP) preprocess-vcf/preprocess_vcf $(DIST_TARGET_DIR)
	$(CP) vcf2multialign/vcf2multialign $(DIST_TARGET_DIR)
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/swift-corelibs-libdispatch/LICENSE $(DIST_TARGET_DIR)/swift-corelibs-libdispatch-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

lib/msa2dag/lib/libMsa2Dag.a:
	$(MAKE) -C lib/msa2dag CXX=$(CXX)

lib/libbio/local.mk: local.mk
	$(CP) local.mk lib/libbio

lib/libbio/src/libbio.a: lib/libbio/local.mk
	$(MAKE) -C lib/libbio

lib/swift-corelibs-libdispatch/build/src/libdispatch.a:
	$(RM) -rf lib/swift-corelibs-libdispatch/build && \
	cd lib/swift-corelibs-libdispatch && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-G Ninja \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DCMAKE_C_FLAGS="$(LIBDISPATCH_CFLAGS)" \
		-DCMAKE_CXX_FLAGS="$(LIBDISPATCH_CXXFLAGS)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(NINJA) -v
