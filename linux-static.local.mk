# Build a mostly static binary.

CC					= gcc-6
CXX					= g++-6
MKDIR				= mkdir
GENGETOPT			= gengetopt
RAGEL				= ragel

BOOST_INCLUDE		= -I/home/tnorri/local/boost-1-60-0-gcc/include
BOOST_LIBS			= -L/home/tnorri/local/boost-1-60-0-gcc/lib -lboost_iostreams
LIBDISPATCH_LIBS	= ../lib/libdispatch/libdispatch-build/src/libdispatch.a ../lib/libpwq/libpwq-build/libpthread_workqueue.a /usr/lib/libkqueue.a -lpthread -static-libstdc++ -static-libgcc
