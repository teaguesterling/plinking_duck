// plinking_pgen_vfs — a thin VFS shim for pgenlib's file opens.
//
// pgenlib opens .pgen/.pgi with raw fopen() and then uses only stdio
// (fseeko/ftello/fread_unlocked/rewind/fclose). The milestone-1 poison enumeration
// proved the op set touching the file is pure stdio, so a real FILE* backed by
// custom callbacks (glibc/musl fopencookie, macOS/BSD funopen) is transparent to
// pgenlib — every read/seek/tell/close works unchanged.
//
// The only pgenlib patch is redirecting its 3 fopen() sites to
// plinking_pgen_fopen(). By default that is a plain fopen (Path L — local/native
// fastpath, byte-identical to unpatched plink2). When the extension registers a
// per-thread VFS opener (Path V), a fname the opener claims is served by a cookie
// FILE* whose reads go through a DuckDB FileHandle (positioned reads) — so remote
// (.pgen over s3://, https://, …) and VFS-resolved paths work with no other
// pgenlib change and no DuckDB types inside pgenlib.
#pragma once

#include <cstdint>
#include <cstdio>

extern "C" {

// A per-thread, context-carrying byte source the extension supplies for Path V.
// `open` returns an opaque handle for `fname` (or nullptr to fall back to fopen);
// the extension owns it (typically a DuckDB FileHandle) and its closures capture
// the ClientContext, so credentials/secrets/settings are handled extension-side.
struct PlinkingPgenVfsOpener {
	void *user;                                                    // opaque extension state
	void *(*open)(const char *fname, void *user);                  // -> handle, or nullptr
	int64_t (*pread)(void *handle, void *buf, int64_t n, uint64_t offset); // bytes, <0 on error
	int64_t (*size)(void *handle);                                 // file size, <0 on error
	void (*close)(void *handle);
};

// pgenlib calls this in place of fopen(fname, "rb"). Path L (default) = real fopen;
// Path V = a cookie FILE* over the registered opener when it claims `fname`.
FILE *plinking_pgen_fopen(const char *fname);

// Register / clear the per-thread VFS opener. The extension sets it (with the
// query's context captured in the callbacks) before initializing a pgen reader on
// a remote/VFS path, and clears it after. Nullptr => Path L (fopen) for all opens.
// The pointer must outlive any pgen open made while it is set.
void plinking_pgen_set_vfs_opener(const PlinkingPgenVfsOpener *opener);

// True iff this platform has the cookie/funopen primitive (Path V available).
int plinking_pgen_vfs_supported(void);

} // extern "C"
