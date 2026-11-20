#include "pgen_header_cache.hpp"
#include "pgen_vfs_opener.hpp"

#include "duckdb/common/file_system.hpp"

#include <mutex>
#include <sys/stat.h>
#include <unordered_map>

namespace duckdb {

std::atomic<int64_t> g_pgen_header_decodes {0};

static shared_ptr<PgenHeader> DecodePgenHeader(ClientContext &context, bool use_vfs, const string &pgen_path,
                                               uint32_t variant_ct, uint32_t sample_ct, const string &what) {
	// No PgenVfsScope here: the caller holds one, and a nested scope's destructor clears the global
	// opener (pgen_vfs_opener.cpp), so the caller's next open (PgrInit) would bypass the VFS.
	auto h = make_shared_ptr<PgenHeader>();
	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;
	plink2::PglErr err = plink2::PgfiInitPhase1(pgen_path.c_str(), nullptr, variant_ct, sample_ct, &header_ctrl,
	                                            &h->pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		throw IOException("%s (phase 1) for '%s': %s", what, pgen_path, errstr_buf);
	}
	if (pgfi_alloc_cacheline_ct > 0) {
		h->alloc.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	}
	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, h->pgfi.raw_variant_ct, &h->max_vrec_width, &h->pgfi,
	                             h->alloc.As<unsigned char>(), &h->pgr_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		throw IOException("%s (phase 2) for '%s': %s", what, pgen_path, errstr_buf);
	}
	// Every reader opens its own handle, so no PgrInit writes to the shared struct.
	if (h->pgfi.shared_ff) {
		fclose(h->pgfi.shared_ff);
		h->pgfi.shared_ff = nullptr;
	}
	h->bytes = pgfi_alloc_cacheline_ct * plink2::kCacheline + sizeof(PgenHeader);
	g_pgen_header_decodes.fetch_add(1, std::memory_order_relaxed);
	return h;
}

//! The cache key: the path and the file's identity. The counts a caller expects are not part of it:
//! a bind that passes UINT32_MAX and a scan that passes the exact counts share one entry, and a hit
//! whose counts differ from the caller's is decoded afresh, uncached, so phase 1 raises pgenlib's own
//! error (AcquirePgenHeader). A stale header gives WRONG records, not lost rows, so the identity is as
//! strict as the file allows:
//! - a path stat() can see (a local file, on any route, including plinking_pgen_io = 'vfs'): device,
//!   inode, size and mtime in nanoseconds. A rewrite that promotes by rename (our writer) gets a new
//!   inode, even within the same second at the same size;
//! - otherwise (a remote path): DuckDB's size, last-modified time and version tag (the ETag, for
//!   httpfs). A file with neither a version tag nor a modification time is not cached.
//! Returns "" when the identity cannot be read: that open is then not cached.
static string PgenHeaderKey(ClientContext &context, const string &pgen_path) {
	string identity;
	struct stat st;
	if (stat(pgen_path.c_str(), &st) == 0) {
		identity = StringUtil::Format("dev=%llu ino=%llu size=%lld mtime=%lld.%09ld", (unsigned long long)st.st_dev,
		                              (unsigned long long)st.st_ino, (long long)st.st_size,
		                              (long long)st.st_mtim.tv_sec, (long)st.st_mtim.tv_nsec);
	} else {
		try {
			auto &fs = FileSystem::GetFileSystem(context);
			auto handle = fs.OpenFile(pgen_path, FileFlags::FILE_FLAGS_READ);
			const string tag = fs.GetVersionTag(*handle);
			const int64_t mtime = fs.GetLastModifiedTime(*handle).value;
			if (tag.empty() && mtime == 0) {
				return "";
			}
			identity = StringUtil::Format("size=%lld mtime=%lld tag=%s", (long long)handle->GetFileSize(),
			                              (long long)mtime, tag);
		} catch (std::exception &) {
			return "";
		}
	}
	return StringUtil::Format("plinking_pgen_header:%s|%s", pgen_path, identity);
}

static bool PgenHeaderCacheEnabled(ClientContext &context) {
	Value v;
	return !context.TryGetCurrentSetting("plinking_pgen_header_cache", v) || v.IsNull() || v.GetValue<bool>();
}

//! Whether a header matches the counts a caller expects (UINT32_MAX = unknown, matches any).
static bool CountsMatch(const PgenHeader &h, uint32_t variant_ct, uint32_t sample_ct) {
	return (variant_ct == UINT32_MAX || variant_ct == h.pgfi.raw_variant_ct) &&
	       (sample_ct == UINT32_MAX || sample_ct == h.pgfi.raw_sample_ct);
}

shared_ptr<PgenHeader> AcquirePgenHeader(ClientContext &context, bool use_vfs, const string &pgen_path,
                                         uint32_t variant_ct, uint32_t sample_ct, const string &what) {
	auto uncached = [&]() {
		return DecodePgenHeader(context, use_vfs, pgen_path, variant_ct, sample_ct, what);
	};
	if (!PgenHeaderCacheEnabled(context) || IsLocalizedPgenTemp(pgen_path)) {
		return uncached();
	}
	const string key = PgenHeaderKey(context, pgen_path);
	if (key.empty()) {
		return uncached();
	}
	auto &cache = ObjectCache::GetObjectCache(context);
	if (auto hit = cache.Get<PgenHeader>(key)) {
		return CountsMatch(*hit, variant_ct, sample_ct) ? hit : uncached();
	}
	// One decode per key: the first caller decodes under the key's lock, later callers wait and share
	// it. Different keys decode in parallel (a full scan opens several cells at once).
	static std::mutex locks_mutex;
	static std::unordered_map<string, shared_ptr<std::mutex>> locks;
	shared_ptr<std::mutex> key_lock;
	{
		std::lock_guard<std::mutex> g(locks_mutex);
		auto &slot = locks[key];
		if (!slot) {
			slot = make_shared_ptr<std::mutex>();
		}
		key_lock = slot;
	}
	// Drop the key's lock from the map on every exit, a throwing decode included, unless another caller
	// still waits on it (the map's reference and ours are the only two then).
	struct EraseLock {
		const string &key;
		~EraseLock() {
			std::lock_guard<std::mutex> g(locks_mutex);
			auto it = locks.find(key);
			if (it != locks.end() && it->second.use_count() <= 2) {
				locks.erase(it);
			}
		}
	} erase_lock {key};
	std::lock_guard<std::mutex> g(*key_lock);
	auto h = cache.Get<PgenHeader>(key);
	if (!h) {
		h = DecodePgenHeader(context, use_vfs, pgen_path, variant_ct, sample_ct, what);
		cache.Put(key, h);
		return h;
	}
	return CountsMatch(*h, variant_ct, sample_ct) ? h : uncached();
}

void InitPgenReaderOn(const PgenHeader &header, const string &pgen_path, plink2::PgenReader &pgr,
                      AlignedBuffer &pgr_alloc, const string &what) {
	if (header.pgr_alloc_cacheline_ct > 0) {
		pgr_alloc.Allocate(header.pgr_alloc_cacheline_ct * plink2::kCacheline);
	} else {
		pgr_alloc.Reset();
	}
	plink2::PreinitPgr(&pgr);
	// PgrInit takes a non-const pointer; with shared_ff closed it only reads the struct (it copies it
	// into pgr.fi and opens pgen_path itself).
	auto err = plink2::PgrInit(pgen_path.c_str(), header.max_vrec_width,
	                           const_cast<plink2::PgenFileInfo *>(&header.pgfi), &pgr, pgr_alloc.As<unsigned char>());
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr ce = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&pgr, &ce);
		throw IOException("%s: PgrInit failed for '%s'", what, pgen_path);
	}
}

} // namespace duckdb
