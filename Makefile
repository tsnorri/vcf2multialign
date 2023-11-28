include local.mk
include common.mk

DEPENDENCIES = lib/libbio/src/libbio.a

# “$() $()” is a literal space.
OS_NAME = $(shell tools/os_name.sh)
VERSION = $(subst $() $(),-,$(shell tools/git_version.sh))
DIST_TARGET_DIR = vcf2multialign-$(VERSION)
DIST_NAME_SUFFIX = $(if $(TARGET_TYPE),-$(TARGET_TYPE),)
DIST_TAR_GZ = vcf2multialign-$(VERSION)-$(OS_NAME)$(DIST_NAME_SUFFIX).tar.gz


.PHONY: all clean-all clean clean-dependencies dependencies

all:	libvcf2multialign/libvcf2multialign.a \
		vcf2multialign/vcf2multialign

clean-all: clean clean-dependencies clean-dist
	$(MAKE) -C tests clean

clean:
	$(MAKE) -C libvcf2multialign clean
	$(MAKE) -C vcf2multialign clean

clean-dependencies: lib/libbio/local.mk
	$(MAKE) -C lib/libbio clean-all

clean-dist:
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

dist: $(DIST_TAR_GZ)

tests: tests/coverage/index.html

libvcf2multialign/libvcf2multialign.a: $(DEPENDENCIES)
	$(MAKE) -C libvcf2multialign

libvcf2multialign/libvcf2multialign.coverage.a: $(DEPENDENCIES)
	$(MAKE) -C libvcf2multialign libvcf2multialign.coverage.a

vcf2multialign/vcf2multialign: $(DEPENDENCIES) libvcf2multialign/libvcf2multialign.a
	$(MAKE) -C vcf2multialign

lib/libbio/vcfcat/vcfcat: $(DEPENDENCIES)
	$(MAKE) -C lib/libbio/vcfcat

$(DIST_TAR_GZ):	vcf2multialign/vcf2multialign
	$(MKDIR) -p $(DIST_TARGET_DIR)
	$(CP) vcf2multialign/vcf2multialign $(DIST_TARGET_DIR)
	$(CP) lib/libbio/vcfcat/vcfcat $(DIST_TARGET_DIR)
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/cereal/LICENSE $(DIST_TARGET_DIR)/cereal-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

lib/libbio/local.mk: local.mk
	$(CP) local.mk lib/libbio

lib/libbio/src/libbio.a: lib/libbio/local.mk
	$(MAKE) -C lib/libbio

lib/libbio/lib/rapidcheck/build/librapidcheck.a:
	$(MAKE) -C lib/libbio lib/rapidcheck/build/librapidcheck.a

tests/coverage/index.html: lib/libbio/lib/rapidcheck/build/librapidcheck.a libvcf2multialign/libvcf2multialign.coverage.a $(DEPENDENCIES)
	$(MAKE) -C tests coverage
