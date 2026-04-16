#pragma once

#include "duckdb.hpp"
#include "duckdb/function/table_function.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

// ---------------------------------------------------------------------------
// SampleInfo — reusable sample metadata for read_pgen / read_pfile (P1-003+)
// ---------------------------------------------------------------------------

//! Sample metadata extracted from .psam or .fam files.
//!
//! At biobank scale (~7M samples), the iid_to_idx map is expensive to build
//! (string hashing + allocations) and most queries do not need it (integer
//! sample filters, or no sample filter at all). So iid_to_idx is **lazy**:
//! call EnsureIidMap() before using it. Callers that do not need IID-based
//! lookups never pay the cost.
//!
//! Similarly, `iids` may be empty when sample_ct is loaded from parquet
//! metadata alone (cheap constant-time call). Callers that need the IID
//! strings (for column names, VARCHAR lookups, etc.) must ensure iids is
//! populated before access.
struct SampleInfo {
	vector<string> iids;                     //!< Individual IDs in file order (possibly lazy)
	vector<string> fids;                     //!< Family IDs (empty if no FID column)
	idx_t sample_ct;                         //!< Total sample count (always populated)
	unordered_map<string, idx_t> iid_to_idx; //!< IID → file-order index (lazy; call EnsureIidMap)

	//! Build iid_to_idx from iids. Validates uniqueness. Safe to call multiple times.
	//! Throws if duplicate IIDs are found. Requires iids to be populated.
	void EnsureIidMap(const string &source_label = "sample file");
};

//! Parse a .psam or .fam file and return sample metadata.
//! Detects format automatically from the first line.
//! Uses DuckDB's VFS for file access (supports S3, HTTP, etc).
SampleInfo LoadSampleInfo(ClientContext &context, const string &path);

// ---------------------------------------------------------------------------
// Psam header parsing — reusable by read_pgen bind phase
// ---------------------------------------------------------------------------

enum class PsamFormat : uint8_t {
	PSAM_FID, //!< Header starts with #FID — FID column present
	PSAM_IID, //!< Header starts with #IID — no FID column
	FAM       //!< No header line — legacy 6-column .fam format
};

//! Result of parsing the first line of a .psam/.fam file.
struct PsamHeaderInfo {
	PsamFormat format;
	vector<string> column_names;      //!< Column names in file order
	vector<LogicalType> column_types; //!< DuckDB types for each column
};

//! Parse the header (or detect .fam format) from already-read file lines.
//! Avoids re-reading the file when the caller already has the lines.
PsamHeaderInfo ParsePsamHeaderFromLines(const vector<string> &lines, const string &path);

//! Parse the header (or detect .fam format) from the given file path.
//! Uses DuckDB's VFS for file access (supports S3, HTTP, etc).
PsamHeaderInfo ParsePsamHeader(ClientContext &context, const string &path);

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

//! Register the read_psam table function with the extension loader.
void RegisterPsamReader(ExtensionLoader &loader);

} // namespace duckdb
