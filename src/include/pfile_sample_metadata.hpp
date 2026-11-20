#pragma once

#include "duckdb.hpp"
#include "psam_reader.hpp"

namespace duckdb {

//! Extended sample info for genotype orient mode: includes all psam columns.
struct PfileSampleMetadata {
	PsamHeaderInfo header;
	//! All data rows from the .psam/.fam, in file order (TEXT path only).
	//! Each inner vector has one field per column. Empty for the parquet path, whose
	//! column values are loaded projection-aware at InitGlobal (PfileGlobalState).
	vector<vector<string>> rows;

	//! Source path for the parquet/non-native path — the projected column data is
	//! (re-)queried from here at InitGlobal, not materialized at bind. Empty for text.
	string source_path;
	//! Sample-sharded cohort over parquet companions: every shard's .psam.parquet, in shard
	//! order (source_paths[0] == source_path). Empty for a single source. See PsamSelectSql.
	vector<string> source_paths;

	//! IID / FID column positions in the header (for the parquet path).
	idx_t iid_col_idx = DConstants::INVALID_INDEX;
	idx_t fid_col_idx = DConstants::INVALID_INDEX;

	//! Index of SEX column in the header (DConstants::INVALID_INDEX if absent)
	idx_t sex_col_idx = DConstants::INVALID_INDEX;
	//! Indices of PAT/MAT columns for special missing-value handling
	vector<idx_t> parent_col_indices;

	//! True when the psam came from a parquet/non-native source: column data is not
	//! in `rows`; the scan reads it from the projection-aware CDC in the global state.
	bool from_parquet = false;
};

} // namespace duckdb
