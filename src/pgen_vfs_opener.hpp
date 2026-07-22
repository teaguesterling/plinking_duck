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
//! always V; localize => not yet implemented (errors).
bool PgenIoUseVfs(ClientContext &context, const string &pgen_path);

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
