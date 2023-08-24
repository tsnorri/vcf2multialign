// This workaround is needed b.c. sysroot_linux-64 from conda-forge (version 2.17)
// contains glibc headers that contain a bug fixed in 2013.
// Details about the fix are available at http://mackyle.github.io/blocksruntime/#glibc
#ifndef UNISTD_FIX_H
#	define UNISTD_FIX_H
#	undef __block
#	include_next <unistd.h>
#	define __block __attribute__((__blocks__(byref)))
#endif
