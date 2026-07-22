// Implementation of the pgenlib VFS shim (see plinking_pgen_vfs.hpp).
//
// Path L (default): plinking_pgen_fopen -> fopen(fname, "rb").
// Path V (opener registered + claims fname): a cookie FILE* (fopencookie on
// glibc/musl, funopen on macOS/BSD) whose read/seek/close drive the extension's
// per-thread opener (positioned reads on a DuckDB FileHandle). glibc stdio buffers
// the cookie stream, so Path V gets readahead for free on top of the handle.
//
// Unsupported platforms (Windows/wasm — already excluded from the build) have no
// cookie primitive: plinking_pgen_vfs_supported() returns 0 and Path V opens fall
// back to fopen (remote simply fails to open, as today), so nothing regresses.

// fopencookie / off64_t / cookie_io_functions_t are GNU extensions.
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <cstdio>

#include "plinking_pgen_vfs.hpp"

#include <cstdlib>
#include <cstring>

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__)
#define PLINKING_PGEN_VFS_FUNOPEN 1
#elif defined(__linux__) || defined(__GLIBC__)
#define PLINKING_PGEN_VFS_COOKIE 1
#endif

namespace {

// Per-thread registered opener (set by the extension for Path V).
thread_local const PlinkingPgenVfsOpener *t_opener = nullptr;

#if defined(PLINKING_PGEN_VFS_COOKIE) || defined(PLINKING_PGEN_VFS_FUNOPEN)
// Cookie state — SELF-CONTAINED: it COPIES the read/close callbacks (it must
// outlive the PlinkingPgenVfsOpener, which the extension may destroy right after
// the open, while the FILE* lives for the whole scan).
struct PgenCookie {
	int64_t (*pread)(void *handle, void *buf, int64_t n, uint64_t offset);
	void (*close)(void *handle);
	void *handle;
	uint64_t offset;
	uint64_t size;
};

PgenCookie *MakeCookie(const PlinkingPgenVfsOpener *ops, void *handle) {
	auto *c = static_cast<PgenCookie *>(std::malloc(sizeof(PgenCookie)));
	if (!c) {
		return nullptr;
	}
	c->pread = ops->pread;
	c->close = ops->close;
	c->handle = handle;
	c->offset = 0;
	int64_t sz = ops->size(handle);
	c->size = (sz < 0) ? 0 : static_cast<uint64_t>(sz);
	return c;
}
#endif

#if defined(PLINKING_PGEN_VFS_COOKIE)
ssize_t CookieRead(void *cookie, char *buf, size_t n) {
	auto *c = static_cast<PgenCookie *>(cookie);
	// Clamp to available (short read at EOF — pgenlib/stdio may request past end).
	uint64_t avail = (c->offset < c->size) ? (c->size - c->offset) : 0;
	if (static_cast<uint64_t>(n) > avail) {
		n = static_cast<size_t>(avail);
	}
	if (n == 0) {
		return 0;
	}
	int64_t got = c->pread(c->handle, buf, static_cast<int64_t>(n), c->offset);
	if (got < 0) {
		return -1;
	}
	c->offset += static_cast<uint64_t>(got);
	return static_cast<ssize_t>(got);
}

int CookieSeek(void *cookie, off64_t *offset, int whence) {
	auto *c = static_cast<PgenCookie *>(cookie);
	int64_t base = (whence == SEEK_SET)   ? 0
	               : (whence == SEEK_CUR) ? static_cast<int64_t>(c->offset)
	                                      : static_cast<int64_t>(c->size);
	int64_t target = base + *offset;
	if (target < 0) {
		return -1;
	}
	c->offset = static_cast<uint64_t>(target);
	*offset = target;
	return 0;
}

int CookieClose(void *cookie) {
	auto *c = static_cast<PgenCookie *>(cookie);
	c->close(c->handle);
	std::free(c);
	return 0;
}

FILE *OpenCookieFile(const PlinkingPgenVfsOpener *ops, void *handle) {
	PgenCookie *c = MakeCookie(ops, handle);
	if (!c) {
		ops->close(handle);
		return nullptr;
	}
	cookie_io_functions_t io = {CookieRead, nullptr, CookieSeek, CookieClose};
	FILE *f = fopencookie(c, "rb", io);
	if (!f) {
		ops->close(handle);
		std::free(c);
	}
	return f;
}
#elif defined(PLINKING_PGEN_VFS_FUNOPEN)
int FunReadFn(void *cookie, char *buf, int n) {
	auto *c = static_cast<PgenCookie *>(cookie);
	uint64_t avail = (c->offset < c->size) ? (c->size - c->offset) : 0;
	if (n < 0) {
		return -1;
	}
	if (static_cast<uint64_t>(n) > avail) {
		n = static_cast<int>(avail);
	}
	if (n == 0) {
		return 0;
	}
	int64_t got = c->pread(c->handle, buf, n, c->offset);
	if (got < 0) {
		return -1;
	}
	c->offset += static_cast<uint64_t>(got);
	return static_cast<int>(got);
}

fpos_t FunSeekFn(void *cookie, fpos_t offset, int whence) {
	auto *c = static_cast<PgenCookie *>(cookie);
	int64_t base = (whence == SEEK_SET)   ? 0
	               : (whence == SEEK_CUR) ? static_cast<int64_t>(c->offset)
	                                      : static_cast<int64_t>(c->size);
	int64_t target = base + static_cast<int64_t>(offset);
	if (target < 0) {
		return -1;
	}
	c->offset = static_cast<uint64_t>(target);
	return static_cast<fpos_t>(target);
}

int FunCloseFn(void *cookie) {
	auto *c = static_cast<PgenCookie *>(cookie);
	c->close(c->handle);
	std::free(c);
	return 0;
}

FILE *OpenCookieFile(const PlinkingPgenVfsOpener *ops, void *handle) {
	PgenCookie *c = MakeCookie(ops, handle);
	if (!c) {
		ops->close(handle);
		return nullptr;
	}
	FILE *f = funopen(c, FunReadFn, nullptr, FunSeekFn, FunCloseFn);
	if (!f) {
		ops->close(handle);
		std::free(c);
	}
	return f;
}
#endif

} // namespace

extern "C" {

FILE *plinking_pgen_fopen(const char *fname) {
#if defined(PLINKING_PGEN_VFS_COOKIE) || defined(PLINKING_PGEN_VFS_FUNOPEN)
	const PlinkingPgenVfsOpener *op = t_opener;
	if (op && op->open) {
		void *handle = op->open(fname, op->user);
		if (handle) {
			return OpenCookieFile(op, handle); // Path V
		}
		// opener declined this fname -> fall through to Path L
	}
#endif
	return std::fopen(fname, "rb"); // Path L (local/native fastpath)
}

void plinking_pgen_set_vfs_opener(const PlinkingPgenVfsOpener *opener) {
	t_opener = opener;
}

int plinking_pgen_vfs_supported(void) {
#if defined(PLINKING_PGEN_VFS_COOKIE) || defined(PLINKING_PGEN_VFS_FUNOPEN)
	return 1;
#else
	return 0;
#endif
}

} // extern "C"
