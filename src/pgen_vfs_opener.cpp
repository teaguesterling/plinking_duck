#include "pgen_vfs_opener.hpp"
#include "plinking_pgen_vfs.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/string_util.hpp"

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

bool PgenIoUseVfs(ClientContext &context, const string &pgen_path) {
	string policy = "auto";
	Value v;
	if (context.TryGetCurrentSetting("plinking_pgen_io", v) && !v.IsNull()) {
		policy = StringUtil::Lower(v.ToString());
	}
	if (policy == "native") {
		return false;
	}
	if (policy == "vfs") {
		return true;
	}
	if (policy == "localize") {
		throw InvalidInputException(
		    "plinking_pgen_io := 'localize' is not yet implemented (use 'auto', 'native', or 'vfs')");
	}
	if (policy != "auto") {
		throw InvalidInputException("unknown plinking_pgen_io := '%s' (expected 'auto', 'native', 'vfs', 'localize')",
		                            policy);
	}
	// auto: route remote/VFS paths through Path V; plain-local uses native fopen.
	return FileSystem::IsRemoteFile(pgen_path);
}

} // namespace duckdb
