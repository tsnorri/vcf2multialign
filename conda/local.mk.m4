GCC_ROOT		= CONDA_PREFIX

CC				= $(GCC_ROOT)/bin/gcc
CXX				= $(GCC_ROOT)/bin/g++

CPPFLAGS		=	--sysroot CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot \
					-isystem =/../../include \
					-DLIBBIO_NO_DISPATCH -DLIBBIO_NO_SAM_READER

LDFLAGS			=	--sysroot CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot \
					-L CONDA_PREFIX/lib \
					-pthread -ldl

BOOST_INCLUDE	=
BOOST_LIBS		= -lboost_iostreams
