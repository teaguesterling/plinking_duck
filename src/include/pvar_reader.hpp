#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Result of parsing a .pvar or .bim file header.
//! Contains schema information and format metadata needed by both
//! read_pvar (P1-001) and read_pgen (P1-003).
struct PvarHeaderInfo {
	//! Column names in normalized output order
	vector<string> column_names;
	//! Column types matching column_names
	vector<LogicalType> column_types;
	//! True if the file is legacy .bim format (no header line)
	bool is_bim;
	//! Number of lines to skip before data begins (comments + header)
	idx_t skip_lines;
};

//! Parse the header of a .pvar or .bim file.
//!
//! Detects format (.pvar vs .bim), extracts column schema, and determines
//! where data lines begin. For .pvar files, skips ## comment lines and
//! parses the #CHROM header. For .bim files, uses the fixed 6-column schema
//! with column order normalization.
//!
//! This function is factored out for reuse by read_pgen (P1-003), which
//! needs variant metadata schema during its bind phase.
PvarHeaderInfo ParsePvarHeader(ClientContext &context, const string &file_path);

//! Register the read_pvar table function with DuckDB.
void RegisterPvarReader(ExtensionLoader &loader);

} // namespace duckdb
