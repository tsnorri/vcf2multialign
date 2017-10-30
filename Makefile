include local.mk
include common.mk

DEPENDENCIES = lib/lemon/build/lemon/libemon.a
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
	$(RM) -r lib/libdispatch/libdispatch-build
	$(RM) -r lib/libpwq/libpwq-build
	$(RM) -r lib/lemon/build/lemon/libemon.a

dependencies: $(DEPENDENCIES)

lib/libdispatch/libdispatch-build/src/libdispatch.a: lib/libpwq/libpwq-build/libpthread_workqueue.a
	rm -rf lib/libdispatch/libdispatch-build && \
	cd lib/libdispatch && \
	mkdir libdispatch-build && \
	cd libdispatch-build && \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	../configure --cc="$(CC)" --c++="$(CXX)" --release -- \
		-DPTHREAD_WORKQUEUE_INCLUDE_DIRS=../../libpwq/include \
		-DPTHREAD_WORKQUEUE_LIBRARIES=../../libpwq/libpwq-build/libpthread_workqueue.a \
		-DBLOCKS_RUNTIME_LIBRARIES=""
	$(MAKE) -C lib/libdispatch/libdispatch-build VERBOSE=1

lib/libpwq/libpwq-build/libpthread_workqueue.a:
	rm -rf lib/libpwq/libpwq-build && \
	cd lib/libpwq && \
	mkdir libpwq-build && \
	cd libpwq-build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	cmake -DSTATIC_WORKQUEUE=ON ..
	$(MAKE) -C lib/libpwq/libpwq-build VERBOSE=1

lib/lemon/build/lemon/libemon.a:
	rm -rf lib/lemon/build && \
	cd lib/lemon && \
	mkdir build && \
	cd build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	cmake ..
	$(MAKE) -C lib/lemon/build VERBOSE=1
	cd lib/lemon/lemon && cp ../build/lemon/config.h ./
