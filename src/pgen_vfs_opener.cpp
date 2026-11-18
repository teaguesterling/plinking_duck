#include "pgen_vfs_opener.hpp"
#include "plinking_pgen_vfs.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/config.hpp"
#include "duckdb/main/client_context_state.hpp"
#include "duckdb/common/unordered_map.hpp"
#include "duckdb_compat.hpp" // CompatRegisteredState

#include <algorithm>
#include <atomic>
#include <cstring>
#include <mutex>
#include <string>

#ifndef _WIN32
#include <unistd.h> // getpid
#else
#include <process.h>
#define getpid _getpid
#endif

namespace duckdb {

// --- Per-query store of the blocks that header opens read (issue #121) -----------
// Each VFS open has its own cookie block cache, so the scan's open re-fetched every block that
// the bind (or InitGlobal) header open had read: over HTTP, a second GET of the same bytes, per
// file, per query, and once more per scan thread. A HEADER open adds each block it reads here;
// every open of the same path takes a stored block from memory. The store lives for one query
// (cleared at QueryEnd), so a file rewritten between two queries is read fresh.
class PgenHeaderBlockStore : public ClientContextState {
public:
	static constexpr const char *KEY = "plinking_duck_pgen_header_blocks";
	//! A bound on the store's memory, per query. A header open past it adds nothing; the scan then
	//! fetches those blocks itself, as before.
	static constexpr idx_t CAP_BYTES = 64ULL << 20; // 64 MiB

	bool TryGet(const string &key, void *buf, int64_t n) {
		std::lock_guard<std::mutex> guard(lock);
		auto it = blocks.find(key);
		if (it == blocks.end()) {
			return false;
		}
		std::memcpy(buf, it->second.data(), static_cast<size_t>(n));
		return true;
	}

	void Put(const string &key, const void *buf, int64_t n) {
		std::lock_guard<std::mutex> guard(lock);
		if (bytes + static_cast<idx_t>(n) > CAP_BYTES || blocks.count(key)) {
			return;
		}
		auto *p = static_cast<const char *>(buf);
		blocks.emplace(key, string(p, p + n));
		bytes += static_cast<idx_t>(n);
	}

	void QueryEnd(ClientContext &context, optional_ptr<ErrorData> error) override {
		std::lock_guard<std::mutex> guard(lock);
		blocks.clear();
		bytes = 0;
	}

private:
	std::mutex lock;
	//! Key: path, NUL, offset, NUL, length (the length differs only at the end of the file).
	unordered_map<string, string> blocks;
	idx_t bytes = 0;
};

// --- The C-ABI opener callbacks: a void* handle is a PgenVfsHandle* ----------------
// user payload = the FileSystem to open through (the context's OpenerFileSystem,
// which injects the ClientContext opener -> httpfs secrets/settings automatically).
namespace {

struct PgenVfsUser {
	FileSystem *fs;
	shared_ptr<PgenHeaderBlockStore> store;
	bool header; // PgenVfsRole::HEADER: add the blocks this open reads to the store
};

//! One open file. It holds its own reference to the store, because a scan reader's FILE* outlives
//! the PgenVfsScope that opened it.
struct PgenVfsHandle {
	unique_ptr<FileHandle> fh;
	string path;
	shared_ptr<PgenHeaderBlockStore> store;
	bool header;
};

void *PgenVfsOpen(const char *fname, void *user) {
	auto *u = static_cast<PgenVfsUser *>(user);
	try {
		auto handle = u->fs->OpenFile(fname, FileFlags::FILE_FLAGS_READ);
		if (!handle) {
			return nullptr;
		}
		auto *h = new PgenVfsHandle();
		h->fh = std::move(handle);
		h->path = fname;
		h->store = u->store;
		h->header = u->header;
		return h; // ownership passes to the cookie (freed in close)
	} catch (...) {
		return nullptr; // decline -> shim falls back to fopen (fails clearly on a URL)
	}
}

int64_t PgenVfsPread(void *handle, void *buf, int64_t n, uint64_t offset) {
	auto *h = static_cast<PgenVfsHandle *>(handle);
	try {
		string key;
		if (h->store) {
			key = h->path;
			key += '\0';
			key += std::to_string(offset);
			key += '\0';
			key += std::to_string(n);
			if (h->store->TryGet(key, buf, n)) {
				return n;
			}
		}
		// The shim clamps n to the file size, so this exact positioned read fits.
		h->fh->Read(buf, static_cast<idx_t>(n), static_cast<idx_t>(offset));
		if (h->store && h->header) {
			h->store->Put(key, buf, n); // only a read that completed: Read throws on a short read
		}
		return n;
	} catch (...) {
		return -1;
	}
}

int64_t PgenVfsSize(void *handle) {
	auto *h = static_cast<PgenVfsHandle *>(handle);
	try {
		return static_cast<int64_t>(h->fh->GetFileSize());
	} catch (...) {
		return -1;
	}
}

void PgenVfsClose(void *handle) {
	delete static_cast<PgenVfsHandle *>(handle);
}

} // namespace

struct PgenVfsScope::State {
	PgenVfsUser user;
	PlinkingPgenVfsOpener opener;
};

PgenVfsScope::PgenVfsScope(ClientContext &context, bool use_vfs, PgenVfsRole role) : active_(use_vfs) {
	if (!active_) {
		return;
	}
	state_ = make_uniq<State>();
	state_->user.fs = &FileSystem::GetFileSystem(context);
	state_->user.store = CompatRegisteredState(context).GetOrCreate<PgenHeaderBlockStore>(PgenHeaderBlockStore::KEY);
	state_->user.header = role == PgenVfsRole::HEADER;
	state_->opener.user = &state_->user;
	state_->opener.open = PgenVfsOpen;
	state_->opener.pread = PgenVfsPread;
	state_->opener.size = PgenVfsSize;
	state_->opener.close = PgenVfsClose;
	plinking_pgen_set_vfs_opener(&state_->opener);
}

PgenVfsScope::~PgenVfsScope() {
	if (active_) {
		plinking_pgen_set_vfs_opener(nullptr);
	}
}

static string PgenIoPolicy(ClientContext &context) {
	string policy = "auto";
	Value v;
	if (context.TryGetCurrentSetting("plinking_pgen_io", v) && !v.IsNull()) {
		policy = StringUtil::Lower(v.ToString());
	}
	if (policy != "auto" && policy != "native" && policy != "vfs" && policy != "localize") {
		throw InvalidInputException("unknown plinking_pgen_io := '%s' (expected 'auto', 'native', 'vfs', 'localize')",
		                            policy);
	}
	return policy;
}

bool PgenIoLocalizes(ClientContext &context) {
	return PgenIoPolicy(context) == "localize";
}

bool PgenIoUseVfs(ClientContext &context, const string &pgen_path) {
	string policy = PgenIoPolicy(context);
	if (policy == "vfs") {
		return true;
	}
	// native, and localize (path already rewritten to a local temp), read natively.
	if (policy == "native" || policy == "localize") {
		return false;
	}
	return PgenPathNeedsVfs(pgen_path);
}

bool PgenPathNeedsVfs(const string &pgen_path) {
	// Route remote/VFS paths through Path V; plain-local uses native fopen. Any other scheme:// path
	// (file://, a protocol file system) also goes through the VFS: fopen cannot open it.
	return FileSystem::IsRemoteFile(pgen_path) || pgen_path.find("://") != string::npos;
}

// --- PgenLocalizeGuard: RAII cleanup of downloaded temp .pgen copies ------------

void PgenLocalizeGuard::Track(FileSystem &fs, string temp_path) {
	fs_ = &fs;
	temp_paths_.push_back(std::move(temp_path));
}

bool PgenLocalizeGuard::Owns(const string &path) const {
	return std::find(temp_paths_.begin(), temp_paths_.end(), path) != temp_paths_.end();
}

void PgenLocalizeGuard::Discard(const string &path) noexcept {
	auto it = std::find(temp_paths_.begin(), temp_paths_.end(), path);
	if (it == temp_paths_.end()) {
		return;
	}
	try {
		if (fs_->FileExists(path)) {
			fs_->RemoveFile(path);
		}
	} catch (...) {
		// best-effort, as in Cleanup.
	}
	temp_paths_.erase(it);
}

void PgenLocalizeGuard::Cleanup() noexcept {
	if (!fs_) {
		return;
	}
	for (auto &p : temp_paths_) {
		try {
			if (fs_->FileExists(p)) {
				fs_->RemoveFile(p);
			}
		} catch (...) {
			// best-effort: a leftover temp is not worth aborting teardown over.
		}
	}
	temp_paths_.clear();
	fs_ = nullptr;
}

PgenLocalizeGuard::~PgenLocalizeGuard() {
	Cleanup();
}

PgenLocalizeGuard::PgenLocalizeGuard(PgenLocalizeGuard &&other) noexcept
    : fs_(other.fs_), temp_paths_(std::move(other.temp_paths_)) {
	other.fs_ = nullptr;
	other.temp_paths_.clear();
}

PgenLocalizeGuard &PgenLocalizeGuard::operator=(PgenLocalizeGuard &&other) noexcept {
	if (this != &other) {
		Cleanup();
		fs_ = other.fs_;
		temp_paths_ = std::move(other.temp_paths_);
		other.fs_ = nullptr;
		other.temp_paths_.clear();
	}
	return *this;
}

// --- localize: download remote (or any) .pgen to a local temp -------------------

//! Directory for localized temps: plinking_localize_dir if set, else DuckDB's
//! configured temporary_directory, else the current directory.
string LocalizeDir(ClientContext &context) {
	Value v;
	if (context.TryGetCurrentSetting("plinking_localize_dir", v) && !v.IsNull()) {
		string dir = v.ToString();
		if (!dir.empty()) {
			return dir;
		}
	}
	auto &config = DBConfig::GetConfig(context);
	if (!config.options.temporary_directory.empty()) {
		return config.options.temporary_directory;
	}
	return ".";
}

namespace {

//! Process-unique temp name (PID + monotonic counter) so concurrent queries
//! localizing the same URL never collide. Skips any name that already exists on
//! disk — defends against a stale temp leaked by a crashed same-PID process (whose
//! counter would also restart at 0), which would otherwise make FILE_CREATE_NEW
//! fail with a confusing error instead of just picking the next name.
string UniqueTempName(FileSystem &fs, const string &dir) {
	static std::atomic<uint64_t> counter {0};
	string path;
	for (int attempt = 0; attempt < 1024; attempt++) {
		uint64_t n = counter.fetch_add(1);
		string name =
		    "plinking_localize_" + std::to_string(static_cast<long>(getpid())) + "_" + std::to_string(n) + ".pgen";
		path = fs.JoinPath(dir, name);
		if (!fs.FileExists(path)) {
			break;
		}
	}
	return path; // if all 1024 existed (absurd), FILE_CREATE_NEW still errors clearly
}

//! Stream-copy src -> dst through `fs` in chunks (no full-file buffer). Sequential
//! reads over the FileHandle (position advances) work over any backing (local or
//! httpfs); the temp is a fresh local file created with FILE_CREATE_NEW.
void LocalizeCopy(FileSystem &fs, const string &src, const string &dst) {
	auto in = fs.OpenFile(src, FileFlags::FILE_FLAGS_READ);
	if (!in) {
		throw IOException("plinking_pgen_io := 'localize': cannot open '%s'", src);
	}
	auto out = fs.OpenFile(dst, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW);
	if (!out) {
		throw IOException("plinking_pgen_io := 'localize': cannot create temp '%s'", dst);
	}
	const idx_t kChunk = 4 * 1024 * 1024; // 4 MiB
	auto buffer = make_unsafe_uniq_array<data_t>(kChunk);
	int64_t got;
	while ((got = in->Read(buffer.get(), kChunk)) > 0) {
		out->Write(buffer.get(), static_cast<idx_t>(got));
	}
	out->Sync();
}

} // namespace

void LocalizePgenIfRequested(ClientContext &context, string &pgen_path, PgenLocalizeGuard &guard) {
	if (PgenIoPolicy(context) != "localize") {
		return;
	}
	LocalizePgen(context, pgen_path, guard, LocalizeDir(context));
}

void LocalizePgen(ClientContext &context, string &pgen_path, PgenLocalizeGuard &guard, const string &dir) {
	auto &fs = FileSystem::GetFileSystem(context);
	// The dir may not exist yet (e.g. DuckDB's default temp dir is created lazily,
	// only when it first spills). FILE_CREATE_NEW won't make parents, so ensure it.
	// Tolerate a concurrent create from another query (benign EEXIST-style race).
	if (!dir.empty() && dir != "." && !fs.DirectoryExists(dir)) {
		try {
			fs.CreateDirectory(dir);
		} catch (...) {
			if (!fs.DirectoryExists(dir)) {
				throw;
			}
		}
	}
	string temp = UniqueTempName(fs, dir);
	string temp_pgi = temp + ".pgi";
	// Track BEFORE the copy so a partial/failed download is still unlinked if the process unwinds
	// past here. A failed copy also removes its temps at once: the path is not rewritten, so a retry
	// (a prepared statement's next EXECUTE) copies again, and each retry would leave one more temp
	// until the bind data dies.
	guard.Track(fs, temp);
	try {
		LocalizeCopy(fs, pgen_path, temp);
		// Split-index fileset: if the source uses an external index, pgenlib will look
		// for the localized .pgen's derived companion ("<temp>.pgi"). We only copied the
		// .pgen, so also copy the source's ".pgi" companion (standard name is
		// "<pgen>.pgi") to that derived name — otherwise the split-index read fails on
		// the missing index. Embedded-index filesets have no .pgi (skipped).
		string src_pgi = pgen_path + ".pgi";
		if (fs.FileExists(src_pgi)) {
			guard.Track(fs, temp_pgi);
			LocalizeCopy(fs, src_pgi, temp_pgi);
		}
	} catch (...) {
		guard.Discard(temp_pgi);
		guard.Discard(temp);
		throw;
	}
	pgen_path = temp; // downstream opens read the native local temp
}

} // namespace duckdb
