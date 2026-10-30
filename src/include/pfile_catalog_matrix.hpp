#pragma once

#include "duckdb/common/common.hpp" // duckdb::string / duckdb::vector aliases
#include "duckdb/common/file_system.hpp"

namespace duckdb {

// ---------------------------------------------------------------------------
// `pfile_catalog` v2 — the matrix (2-D) model (spec: docs/planning/P13-001-pfile-catalog-matrix §3/§4).
//
// A grid of cells: a cell is one variant-block x one sample-block = one `.pgen`. Companions are
// shared ALONG each axis (§3.1): every cell in a variant-block shares that block's `.pvar`; every
// cell in a sample-block shares that block's `.psam`.
//
// Phase 1 grammar (one `enum` key per axis; `ordinal_interval`/`value_bin` and holes come later):
//   { "pfile_catalog": 2, "dataset": "...",
//     "axes": [ {"name":"variants","keys":[{"name":"CHROM","type":"enum",
//                  "values":[{"value":"1","start":0,"stop":4}, ...]}]},
//               {"name":"samples", "keys":[{"name":"GROUP", ...}]} ],
//     "template": {"pfile":"data/CHROM={CHROM}/GROUP={GROUP}/cell",
//                  "pvar":"ref/CHROM={CHROM}/block.pvar", "psam":"ref/GROUP={GROUP}/block.psam"} }
// `{KEY}` in a template is replaced by the key's value; `pvar` may reference only variants-axis
// keys, `psam` only samples-axis keys (§5 invariant, generalized). Declared spans are half-open
// `[start, stop)` and must tile the axis (§6 boundary law).
// ---------------------------------------------------------------------------

//! One value of an `enum` key and the half-open ordinal span of the block it labels.
struct AxisValue {
	string value;
	int64_t start = 0;
	int64_t stop = 0;
};

//! A key on an axis. Phase 1: `type == "enum"` only.
struct AxisKey {
	string name;
	string type;
	vector<AxisValue> values; //! ordered
};

//! An axis (`variants` | `samples`) and its (phase 1: single) key.
struct AxisDef {
	string name;
	AxisKey key;
};

//! One grid cell, resolved. `keys` is (key name, value) for the variants key then the samples key.
struct MatrixCell {
	vector<std::pair<string, string>> keys;
	string prefix; //! `.pgen` path without extension (absolute, normalized)
	string pvar;   //! the variant-block's shared `.pvar`
	string psam;   //! the sample-block's shared `.psam`
	int64_t variant_start = 0;
	int64_t variant_stop = 0;
	int64_t sample_start = 0;
	int64_t sample_stop = 0;
	int64_t source_ordinal = 0;
};

struct MatrixCatalog {
	int64_t version = 2;
	string dataset;
	string base_dir;
	AxisDef variants;
	AxisDef samples;
	string pfile_template; //! unresolved (relative to base_dir)
	string pvar_template;
	string psam_template;
};

//! True when `content` is a JSON object whose integer `pfile_catalog` is 2 (the version router).
bool IsMatrixCatalog(const string &content);

//! Parse + validate a v2 matrix catalog. Throws "read_pfile_catalog: catalog parse error: ...".
MatrixCatalog ParseMatrixCatalog(FileSystem &fs, const string &path, const string &content);

//! Dense Cartesian product of the two axes' blocks, variants-major; resolves each cell's `.pgen` from
//! the template and its shared `.pvar`/`.psam`. Opens no file.
vector<MatrixCell> BuildCells(const MatrixCatalog &cat, FileSystem &fs);

} // namespace duckdb
