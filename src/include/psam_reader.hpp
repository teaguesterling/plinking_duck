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
//! Used by read_psam for its own output and by read_pgen to look up sample
//! information without re-parsing the file.
struct SampleInfo {
	vector<string> iids;                       //!< Individual IDs in file order
	vector<string> fids;                       //!< Family IDs (empty if no FID column)
	idx_t sample_ct;                           //!< Total sample count
	unordered_map<string, idx_t> iid_to_idx;   //!< IID → file-order index
};

//! Parse a .psam or .fam file and return sample metadata.
//! Detects format automatically from the first line.
SampleInfo LoadSampleInfo(const string &path);

// ---------------------------------------------------------------------------
// Psam header parsing — reusable by read_pgen bind phase
// ---------------------------------------------------------------------------

enum class PsamFormat : uint8_t {
	PSAM_FID,   //!< Header starts with #FID — FID column present
	PSAM_IID,   //!< Header starts with #IID — no FID column
	FAM         //!< No header line — legacy 6-column .fam format
};

//! Result of parsing the first line of a .psam/.fam file.
struct PsamHeaderInfo {
	PsamFormat format;
	vector<string> column_names;       //!< Column names in file order
	vector<LogicalType> column_types;  //!< DuckDB types for each column
};

//! Parse the header (or detect .fam format) from the given file path.
PsamHeaderInfo ParsePsamHeader(const string &path);

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

//! Register the read_psam table function with the extension loader.
void RegisterPsamReader(ExtensionLoader &loader);

} // namespace duckdb
