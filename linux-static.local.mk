# Build a mostly static binary.

LLVM_ROOT				= /usr/lib/llvm-6.0
CLANG_INCLUDE_DIR		= $(LLVM_ROOT)/lib/clang/6.0.0/include

CC						= clang-6.0
CXX						= clang++-6.0
LIBDISPATCH_CFLAGS		= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)
LIBDISPATCH_CXXFLAGS	= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)

CFLAGS					= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)
CXXFLAGS				= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)

BOOST_LIBS				= -L$(BOOST_ROOT)/lib -lboost_iostreams
LDFLAGS					= -static-libstdc++ -static-libgcc

# Used in .tar.gz name.
TARGET_TYPE				= static
