include local.mk
include common.mk

DEPENDENCIES = lib/libbio/src/libbio.a
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

all:	libvcf2multialign/libvcf2multialign.a \
		preprocess-vcf/preprocess_vcf \
		create-variant-graph/create_variant_graph \
		variant-graph-to-founders/variant_graph_to_founders \
		variant-graph-to-gv/variant_graph_to_gv

clean-all: clean clean-dependencies clean-dist

clean:
	$(MAKE) -C libvcf2multialign clean
	$(MAKE) -C preprocess-vcf clean
	$(MAKE) -C create-variant-graph clean
	$(MAKE) -C variant-graph-to-founders clean
	$(MAKE) -C variant-graph-to-gv clean
	$(MAKE) -C inspect-variant-graph clean

clean-dependencies: lib/libbio/local.mk
	$(MAKE) -C lib/libbio clean-all
	$(RM) -r lib/swift-corelibs-libdispatch/build

clean-dist:
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

dist: $(DIST_TAR_GZ)

libvcf2multialign/libvcf2multialign.a: $(DEPENDENCIES)
	$(MAKE) -C libvcf2multialign

create-variant-graph/create_variant_graph: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C create-variant-graph

inspect-variant-graph/inspect_variant_graph: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C inspect-variant-graph

preprocess-vcf/preprocess_vcf: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C preprocess-vcf

variant-graph-to-founders/variant_graph_to_founders: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C variant-graph-to-founders

variant-graph-to-gv/variant_graph_to_gv: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C variant-graph-to-gv

$(DIST_TAR_GZ):	preprocess-vcf/preprocess_vcf \
				create-variant-graph/create_variant_graph \
				inspect-variant-graph/inspect_variant_graph \
				variant-graph-to-founders/variant_graph_to_founders \
				variant-graph-to-gv/variant_graph_to_gv
	$(MKDIR) -p $(DIST_TARGET_DIR)
	$(CP) preprocess-vcf/preprocess_vcf $(DIST_TARGET_DIR)
	$(CP) create-variant-graph/create_variant_graph $(DIST_TARGET_DIR)
	$(CP) inspect-variant-graph/inspect_variant_graph $(DIST_TARGET_DIR)
	$(CP) variant-graph-to-founders/variant_graph_to_founders $(DIST_TARGET_DIR)
	$(CP) variant-graph-to-gv/variant_graph_to_gv $(DIST_TARGET_DIR)
	$(CP) vcf2multialign/vcf2multialign $(DIST_TARGET_DIR)
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/swift-corelibs-libdispatch/LICENSE $(DIST_TARGET_DIR)/swift-corelibs-libdispatch-license.txt
	$(CP) lib/cereal/LICENSE $(DIST_TARGET_DIR)/cereal-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

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
