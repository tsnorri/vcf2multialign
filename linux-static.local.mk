# Build a mostly static binary.

CC					= gcc-6
CXX					= g++-6
MKDIR				= mkdir
GENGETOPT			= gengetopt
RAGEL				= ragel
DOT					= dot

BOOST_INCLUDE		= -I/home/tnorri/local/boost-1-63-0/include
BOOST_LIBS			= -L/home/tnorri/local/boost-1-63-0/lib -lboost_iostreams
LIBDISPATCH_LIBS	= ../lib/libdispatch/libdispatch-build/src/libdispatch.a ../lib/libpwq/libpwq-build/libpthread_workqueue.a /usr/lib/libkqueue.a -lpthread -static-libstdc++ -static-libgcc
