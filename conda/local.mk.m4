PREFIX		= CONDA_PREFIX
CC			= $(PREFIX)/bin/clang
CXX			= $(PREFIX)/bin/clang++

#LIBBSD_LIB	= $(PREFIX)/x86_64-conda-linux-gnu/sysroot/usr/lib64/libbsd.a
LIBBSD_LIB	= 

INCLUDES	= --sysroot $(PREFIX)/x86_64-conda-linux-gnu/sysroot -isystem $(PREFIX)/x86_64-conda-linux-gnu/include/c++/13.2.0 -isystem $(PREFIX)/x86_64-conda-linux-gnu/include/c++/13.2.0/x86_64-conda-linux-gnu -isystem $(PREFIX)/x86_64-conda-linux-gnu/include -isystem $(PREFIX)/include # -I'${SRC_DIR}/conda'"
LIBRARIES	= -L$(PREFIX)/x86_64-conda-linux-gnu/lib -L$(PREFIX)/lib/gcc/x86_64-conda-linux-gnu/13.2.0 -L$(PREFIX)/lib -pthread -lboost_iostreams $(LIBBSD_LIB) -lz -ldl

#INCLUDES	= --sysroot $(PREFIX) -isystem $(PREFIX)/x86_64-conda-linux-gnu/include/c++/13.2.0 -isystem $(PREFIX)/x86_64-conda-linux-gnu/include/c++/13.2.0/x86_64-conda-linux-gnu -isystem $(PREFIX)/x86_64-conda-linux-gnu/sysroot/usr/include -isystem $(PREFIX)/x86_64-conda-linux-gnu/include
#LIBRARIES	= -L$(PREFIX)/x86_64-conda-linux-gnu/lib -L$(PREFIX)/lib/gcc/x86_64-conda-linux-gnu/13.2.0 -L$(PREFIX)/lib

CPPFLAGS	= $(INCLUDES)
CFLAGS		= -fblocks -pthread -Wno-unknown-warning-option -Wno-unused-command-line-argument
CXXFLAGS	= -fblocks -pthread -Wno-unknown-warning-option -Wno-unused-command-line-argument
LDFLAGS		= --sysroot $(PREFIX)/x86_64-conda-linux-gnu/sysroot $(LIBRARIES)

LIBDISPATCH_CFLAGS		= $(CPPFLAGS) -Wno-unused-command-line-argument $(LIBRARIES)
LIBDISPATCH_CXXFLAGS	= $(CPPFLAGS) -Wno-unused-command-line-argument $(LIBRARIES)
BOOST_ROOT				= $(PREFIX)
