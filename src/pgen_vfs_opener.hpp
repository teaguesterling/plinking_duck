// DuckDB-side glue for the pgenlib VFS shim (src/plinking_pgen_vfs.*).
//
// PgenVfsScope registers a per-thread opener so that, while it is active, pgenlib's
// fopen (via plinking_pgen_fopen) returns a cookie FILE* backed by a DuckDB
// FileHandle opened with the query's context (credentials/secrets/settings). This
// is where the only DuckDB dependency lives — the shim itself is DuckDB-free.
#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Resolve the `plinking_pgen_io` policy against a path: should this pgen open go
//! through DuckDB's VFS (Path V, cookie over FileHandle) vs a native local fopen
//! (Path L)? auto (default) => remote paths use V; native => always L; vfs =>
//! always V; localize => L (the path has already been rewritten to a local temp by
//! LocalizePgenIfRequested). Pure: no I/O, safe to call more than once.
bool PgenIoUseVfs(ClientContext &context, const string &pgen_path);

//! RAII owner of localized (downloaded-to-local-temp) .pgen copies. Lives as a
//! bind_data member so the temps survive the whole query and are unlinked when the
//! bind data is destroyed (per-query temp lifecycle). Move-only: two owners must
//! never unlink the same temp. Multi-source readers (sharded read_pfile) track one
//! temp per localized source.
class PgenLocalizeGuard {
public:
	PgenLocalizeGuard() = default;
	~PgenLocalizeGuard();
	PgenLocalizeGuard(PgenLocalizeGuard &&other) noexcept;
	PgenLocalizeGuard &operator=(PgenLocalizeGuard &&other) noexcept;
	PgenLocalizeGuard(const PgenLocalizeGuard &) = delete;
	PgenLocalizeGuard &operator=(const PgenLocalizeGuard &) = delete;

	//! Record a temp path (with the FileSystem that can remove it) to unlink at
	//! destruction. Call BEFORE writing the temp, so a failed/partial copy is still
	//! cleaned up.
	void Track(FileSystem &fs, string temp_path);

private:
	void Cleanup() noexcept;
	FileSystem *fs_ = nullptr;
	vector<string> temp_paths_;
};

//! If `plinking_pgen_io == 'localize'`, download `pgen_path` to a local temp (in
//! `plinking_localize_dir`, else DuckDB's temp dir), rewrite `pgen_path` in place to
//! the temp, and register the temp with `guard` for query-end cleanup. Reads through
//! the context FileSystem so httpfs secrets/settings apply. No-op for every other
//! policy. Call ONCE per source at bind, before any pgen open on that path.
void LocalizePgenIfRequested(ClientContext &context, string &pgen_path, PgenLocalizeGuard &guard);

//! RAII: while in scope AND use_vfs, pgen opens on THIS thread are served from the
//! VFS. A no-op (native fopen) when use_vfs is false. Register/clear is per-thread,
//! so construct one around each pgenlib open region (bind, per-thread reader init).
class PgenVfsScope {
public:
	PgenVfsScope(ClientContext &context, bool use_vfs);
	~PgenVfsScope();
	PgenVfsScope(const PgenVfsScope &) = delete;
	PgenVfsScope &operator=(const PgenVfsScope &) = delete;

private:
	bool active_;
	struct State;
	unique_ptr<State> state_;
};

} // namespace duckdb
