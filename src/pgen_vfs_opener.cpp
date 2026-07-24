#include "pgen_vfs_opener.hpp"
#include "plinking_pgen_vfs.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/config.hpp"

#include <atomic>

#ifndef _WIN32
#include <unistd.h> // getpid
#else
#include <process.h>
#define getpid _getpid
#endif

namespace duckdb {

// --- The C-ABI opener callbacks: a void* handle is a FileHandle* ---------------
// user payload = the FileSystem to open through (the context's OpenerFileSystem,
// which injects the ClientContext opener -> httpfs secrets/settings automatically).
namespace {

struct PgenVfsUser {
	FileSystem *fs;
};

void *PgenVfsOpen(const char *fname, void *user) {
	auto *u = static_cast<PgenVfsUser *>(user);
	try {
		auto handle = u->fs->OpenFile(fname, FileFlags::FILE_FLAGS_READ);
		if (!handle) {
			return nullptr;
		}
		return handle.release(); // ownership passes to the cookie (freed in close)
	} catch (...) {
		return nullptr; // decline -> shim falls back to fopen (fails clearly on a URL)
	}
}

int64_t PgenVfsPread(void *handle, void *buf, int64_t n, uint64_t offset) {
	auto *fh = static_cast<FileHandle *>(handle);
	try {
		// The shim clamps n to the file size, so this exact positioned read fits.
		fh->Read(buf, static_cast<idx_t>(n), static_cast<idx_t>(offset));
		return n;
	} catch (...) {
		return -1;
	}
}

int64_t PgenVfsSize(void *handle) {
	auto *fh = static_cast<FileHandle *>(handle);
	try {
		return static_cast<int64_t>(fh->GetFileSize());
	} catch (...) {
		return -1;
	}
}

void PgenVfsClose(void *handle) {
	delete static_cast<FileHandle *>(handle);
}

} // namespace

struct PgenVfsScope::State {
	PgenVfsUser user;
	PlinkingPgenVfsOpener opener;
};

PgenVfsScope::PgenVfsScope(ClientContext &context, bool use_vfs) : active_(use_vfs) {
	if (!active_) {
		return;
	}
	state_ = make_uniq<State>();
	state_->user.fs = &FileSystem::GetFileSystem(context);
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

bool PgenIoUseVfs(ClientContext &context, const string &pgen_path) {
	string policy = PgenIoPolicy(context);
	if (policy == "vfs") {
		return true;
	}
	// native, and localize (path already rewritten to a local temp), read natively.
	if (policy == "native" || policy == "localize") {
		return false;
	}
	// auto: route remote/VFS paths through Path V; plain-local uses native fopen.
	return FileSystem::IsRemoteFile(pgen_path);
}

// --- PgenLocalizeGuard: RAII cleanup of downloaded temp .pgen copies ------------

void PgenLocalizeGuard::Track(FileSystem &fs, string temp_path) {
	fs_ = &fs;
	temp_paths_.push_back(std::move(temp_path));
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

namespace {

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

//! Process-unique temp name (PID + monotonic counter) so concurrent queries
//! localizing the same URL never collide.
string UniqueTempName(FileSystem &fs, const string &dir) {
	static std::atomic<uint64_t> counter {0};
	uint64_t n = counter.fetch_add(1);
	string name = "plinking_localize_" + std::to_string(static_cast<long>(getpid())) + "_" + std::to_string(n) + ".pgen";
	return fs.JoinPath(dir, name);
}

//! Stream-copy src -> dst through `fs` in chunks (no full-file buffer). Positioned
//! reads keep this correct over any FileHandle (local or httpfs).
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
	auto &fs = FileSystem::GetFileSystem(context);
	string dir = LocalizeDir(context);
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
	// Track BEFORE the copy so a partial/failed download is still unlinked.
	guard.Track(fs, temp);
	LocalizeCopy(fs, pgen_path, temp);
	pgen_path = temp; // downstream opens read the native local temp
}

} // namespace duckdb
