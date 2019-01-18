# Build a mostly static binary.

CC					= gcc-7
CXX					= g++-7

BOOST_ROOT			= /home/tnorri/local/boost-1-64-0-g++-7.1
BOOST_INCLUDE		= -I$(BOOST_ROOT)/include
BOOST_LIBS			= -L$(BOOST_ROOT)/lib -lboost_iostreams
LIBDISPATCH_LIBS	= ../lib/libdispatch/libdispatch-build/src/libdispatch.a ../lib/libpwq/libpwq-build/libpthread_workqueue.a /usr/lib/libkqueue.a -lpthread -static-libstdc++ -static-libgcc
