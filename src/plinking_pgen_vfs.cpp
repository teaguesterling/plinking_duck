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

// Read-ahead block cache. pgenlib reads in ~8 KiB stdio chunks; a remote FileHandle
// (httpfs) over-fetches a large window PER read and does not reuse it, so a naive
// 1:1 mapping fetched ~20x the bytes pgenlib actually needs. Instead we fetch
// aligned BLOCK-sized chunks once and serve pgenlib's small reads (and pgen's
// seek-back ldbase re-reads) from a tiny LRU — collapsing the over-fetch to ~1x
// while staying below cache_httpfs (which, if loaded, caches these block fetches).
static const uint64_t kPgenVfsBlock = 262144; // 256 KiB
static const int kPgenVfsNBlocks = 8;         // 2 MiB read-ahead cache per open (per reader/thread)

// Cookie state — SELF-CONTAINED: it COPIES the read/close callbacks (it must
// outlive the PlinkingPgenVfsOpener, which the extension may destroy right after
// the open, while the FILE* lives for the whole scan).
struct PgenCookie {
	int64_t (*pread)(void *handle, void *buf, int64_t n, uint64_t offset);
	void (*close)(void *handle);
	void *handle;
	uint64_t offset;
	uint64_t size;
	unsigned char *cache;               // kPgenVfsNBlocks * kPgenVfsBlock
	uint64_t blk_off[kPgenVfsNBlocks];  // aligned block start; UINT64_MAX = empty
	uint32_t blk_len[kPgenVfsNBlocks];  // valid bytes in the block
	uint32_t blk_tick[kPgenVfsNBlocks]; // LRU stamp
	uint32_t tick;
};

PgenCookie *MakeCookie(const PlinkingPgenVfsOpener *ops, void *handle) {
	auto *c = static_cast<PgenCookie *>(std::malloc(sizeof(PgenCookie)));
	if (!c) {
		return nullptr;
	}
	c->cache = static_cast<unsigned char *>(std::malloc(static_cast<size_t>(kPgenVfsNBlocks) * kPgenVfsBlock));
	if (!c->cache) {
		std::free(c);
		return nullptr;
	}
	c->pread = ops->pread;
	c->close = ops->close;
	c->handle = handle;
	c->offset = 0;
	int64_t sz = ops->size(handle);
	c->size = (sz < 0) ? 0 : static_cast<uint64_t>(sz);
	for (int i = 0; i < kPgenVfsNBlocks; i++) {
		c->blk_off[i] = UINT64_MAX;
		c->blk_len[i] = 0;
		c->blk_tick[i] = 0;
	}
	c->tick = 0;
	return c;
}

// Serve up to one block's worth from `c->offset`; returns bytes served (0 at EOF,
// <0 on error). A read spanning a block boundary is clamped — stdio re-calls for
// the rest. Advances c->offset by the amount served.
int64_t CachedRead(PgenCookie *c, void *buf, size_t n) {
	if (c->offset >= c->size) {
		return 0;
	}
	uint64_t avail_total = c->size - c->offset;
	if (static_cast<uint64_t>(n) > avail_total) {
		n = static_cast<size_t>(avail_total);
	}
	if (n == 0) {
		return 0;
	}
	uint64_t blk = (c->offset / kPgenVfsBlock) * kPgenVfsBlock;
	size_t in_off = static_cast<size_t>(c->offset - blk);
	size_t max_in_block = static_cast<size_t>(kPgenVfsBlock) - in_off;
	if (n > max_in_block) {
		n = max_in_block;
	}
	int slot = -1;
	for (int i = 0; i < kPgenVfsNBlocks; i++) {
		if (c->blk_off[i] == blk) {
			slot = i;
			break;
		}
	}
	if (slot < 0) {
		slot = 0; // evict LRU
		for (int i = 1; i < kPgenVfsNBlocks; i++) {
			if (c->blk_tick[i] < c->blk_tick[slot]) {
				slot = i;
			}
		}
		uint64_t want = kPgenVfsBlock;
		if (blk + want > c->size) {
			want = c->size - blk;
		}
		int64_t got =
		    c->pread(c->handle, c->cache + static_cast<size_t>(slot) * kPgenVfsBlock, static_cast<int64_t>(want), blk);
		if (got < 0) {
			c->blk_off[slot] = UINT64_MAX;
			return -1;
		}
		c->blk_off[slot] = blk;
		c->blk_len[slot] = static_cast<uint32_t>(got);
	}
	if (in_off >= c->blk_len[slot]) {
		return 0; // requested past the block's valid bytes (EOF)
	}
	size_t served = n;
	if (served > c->blk_len[slot] - in_off) {
		served = c->blk_len[slot] - in_off;
	}
	std::memcpy(buf, c->cache + static_cast<size_t>(slot) * kPgenVfsBlock + in_off, served);
	c->blk_tick[slot] = ++c->tick;
	c->offset += served;
	return static_cast<int64_t>(served);
}
#endif

#if defined(PLINKING_PGEN_VFS_COOKIE)
ssize_t CookieRead(void *cookie, char *buf, size_t n) {
	return static_cast<ssize_t>(CachedRead(static_cast<PgenCookie *>(cookie), buf, n));
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
	std::free(c->cache);
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
		std::free(c->cache);
		std::free(c);
	}
	return f;
}
#elif defined(PLINKING_PGEN_VFS_FUNOPEN)
int FunReadFn(void *cookie, char *buf, int n) {
	if (n < 0) {
		return -1;
	}
	return static_cast<int>(CachedRead(static_cast<PgenCookie *>(cookie), buf, static_cast<size_t>(n)));
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
	std::free(c->cache);
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
		std::free(c->cache);
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
