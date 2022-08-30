# Build a static binary.

LLVM_ROOT				= /usr/lib/llvm-15
CLANG_INCLUDE_DIR		= $(LLVM_ROOT)/lib/clang/15/include

CC							= clang-15
CXX							= clang++-15
LIBDISPATCH_CFLAGS			= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) 
LIBDISPATCH_CXXFLAGS		= -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) -stdlib=libc++

CFLAGS					= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR)
CXXFLAGS				= -fblocks -U__STDC_HOSTED__ -isystem $(CLANG_INCLUDE_DIR) -stdlib=libc++

BOOST_ROOT				= /home/tnorri/local/boost-1.80.0-clang-15
BOOST_LIBS				= -L$(BOOST_ROOT)/lib -lboost_iostreams
LDFLAGS					= -static -static-libgcc -stdlib=libc++ -lc++ -lpthread -lbsd -lz -ldl

# Used in .tar.gz name.
TARGET_TYPE				= static
