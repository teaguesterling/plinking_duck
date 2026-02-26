// musl_compat.h — rawmemchr shim for non-glibc systems (e.g., musl libc)
//
// plink-ng's plink2_string.h calls rawmemchr() directly when _GNU_SOURCE is
// defined, expecting glibc to provide it.  musl defines _GNU_SOURCE features
// but does not provide rawmemchr (a GNU extension, not POSIX).
//
// This header is force-included (-include) into pgenlib compilation units
// when rawmemchr is not available in the system libc.

#ifndef PLINK_MUSL_COMPAT_H
#define PLINK_MUSL_COMPAT_H

#include <string.h>

#ifdef __cplusplus
extern "C"
#endif
    inline void *
    rawmemchr(const void *s, int c) {
	const unsigned char *p = (const unsigned char *)s;
	while (*p != (unsigned char)c)
		++p;
	return (void *)p;
}

#endif // PLINK_MUSL_COMPAT_H
