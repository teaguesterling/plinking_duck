#pragma once

#include "duckdb/common/common.hpp" // duckdb::string / duckdb::vector aliases
#include "duckdb/common/file_system.hpp"

#include "yyjson.hpp"

namespace duckdb {

// ---------------------------------------------------------------------------
// Helpers shared by the two `pfile_catalog` parsers: the 1-D reader (src/pfile_catalog_reader.cpp,
// which defines the path helpers) and the v2 matrix parser (src/pfile_catalog_matrix.cpp).
// ---------------------------------------------------------------------------

//! RAII for a parsed yyjson document, so every throw path frees it.
struct YyDocGuard {
	duckdb_yyjson::yyjson_doc *doc = nullptr;
	~YyDocGuard() {
		if (doc) {
			duckdb_yyjson::yyjson_doc_free(doc);
		}
	}
};

//! Read a JSON integer into int64; false when absent, not an integer, or out of int64 range.
bool JsonGetInt64(duckdb_yyjson::yyjson_val *v, int64_t &out);

//! Collapse `.`/`..`/duplicate separators in a POSIX-style path (os.path.normpath semantics).
string NormPath(const string &path);
//! Absolute-ize `path` against the process CWD (os.path.abspath), then normalize.
string AbsPath(FileSystem &fs, const string &path);
//! Directory portion of a path. "." when the path has no separator.
string DirName(const string &path);
//! Resolve a catalog path relative to `base_dir` (§7.2); result is absolute + normalized.
string ResolvePath(FileSystem &fs, const string &base_dir, const string &path, bool path_is_relative);

} // namespace duckdb
