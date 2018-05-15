include local.mk
include common.mk

DEPENDENCIES = lib/libbio/src/libbio.a
ifeq ($(shell uname -s),Linux)
	DEPENDENCIES    +=  lib/libdispatch/libdispatch-build/src/libdispatch.a
	DEPENDENCIES    +=  lib/libpwq/libpwq-build/libpthread_workqueue.a
endif


.PHONY: all clean-all clean clean-dependencies dependencies

all: dependencies
	$(MAKE) -C src

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C src clean

clean-dependencies:
	$(MAKE) -C lib/libbio clean-all
	$(RM) -r lib/libdispatch/libdispatch-build
	$(RM) -r lib/libpwq/libpwq-build

dependencies: $(DEPENDENCIES)

lib/libbio/src/libbio.a:
	$(CP) local.mk lib/libbio
	$(MAKE) -C lib/libbio

lib/libdispatch/libdispatch-build/src/libdispatch.a: lib/libpwq/libpwq-build/libpthread_workqueue.a
	$(RM) -rf lib/libdispatch/libdispatch-build && \
	cd lib/libdispatch && \
	$(MKDIR) libdispatch-build && \
	cd libdispatch-build && \
	../configure --cc="$(CC)" --c++="$(CXX)" --release -- \
		-DPTHREAD_WORKQUEUE_INCLUDE_DIRS=../../libpwq/include \
		-DPTHREAD_WORKQUEUE_LIBRARIES=../../libpwq/libpwq-build/libpthread_workqueue.a \
		-DBLOCKS_RUNTIME_LIBRARIES=""
	$(MAKE) -C lib/libdispatch/libdispatch-build VERBOSE=1


lib/libpwq/libpwq-build/libpthread_workqueue.a:
	$(RM) -rf lib/libpwq/libpwq-build && \
	cd lib/libpwq && \
	$(MKDIR) libpwq-build && \
	cd libpwq-build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	$(CMAKE) -DSTATIC_WORKQUEUE=ON ..
	$(MAKE) -C lib/libpwq/libpwq-build VERBOSE=1
