GCC_ROOT = CONDA_PREFIX
CLANG_ROOT = CONDA_PREFIX


CC	= $(CLANG_ROOT)/bin/clang
CXX	= $(CLANG_ROOT)/bin/clang++

CPPFLAGS	=	-nostdinc \
				-isystem CONDA_PREFIX/lib/clang/16/include \
				-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include \
				-isystem CONDA_PREFIX/include

#				-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0 \
#				-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0/x86_64-conda-linux-gnu \

CFLAGS			= -fblocks
CXXFLAGS		= -fblocks -nostdinc++
LDFLAGS			= -L CONDA_PREFIX/lib -pthread -lz -lbz2 -lm -lc -ldl -fblocks

SYSTEM_CXXFLAGS	=

LIBDISPATCH_CFLAGS      = 
LIBDISPATCH_CXXFLAGS    = 
LIBDISPATCH_LDFLAGS     = 

BOOST_INCLUDE	=
BOOST_LIBS		= -lboost_iostreams
