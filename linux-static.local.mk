# Build a static binary.

CC						= gcc-12
CXX						= g++-12
GCOV					= gcov-12

CPPFLAGS				= -DBOOST_STACKTRACE_USE_BACKTRACE -DBOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE="</usr/lib/gcc/x86_64-linux-gnu/12/include/backtrace.h>" -DLIBBIO_NO_DISPATCH -DLIBBIO_NO_SAM_READER

BOOST_ROOT				= /usr
BOOST_LIBS				= -L$(BOOST_ROOT)/lib -lboost_iostreams
BOOST_INCLUDE			=
LDFLAGS					= -static -static-libgcc -pthread -lboost_stacktrace_backtrace -lbacktrace -lbsd -lz -ldl

# Used in .tar.gz name.
TARGET_TYPE				= static
