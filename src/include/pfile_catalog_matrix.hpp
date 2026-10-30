#pragma once

#include "duckdb/common/common.hpp" // duckdb::string / duckdb::vector aliases
#include "duckdb/common/file_system.hpp"

#include <limits>

namespace duckdb {

// ---------------------------------------------------------------------------
// `pfile_catalog` v2 — the matrix (2-D) model (spec: docs/planning/P13-001-pfile-catalog-matrix §3/§4).
//
// A grid of cells: a cell is one variant-block x one sample-block = one `.pgen`. Companions are
// shared ALONG each axis (§3.1): every cell in a variant-block shares that block's `.pvar`; every
// cell in a sample-block shares that block's `.psam`.
//
// Grammar — an axis declares an ORDERED list of typed keys (§3.2) and the BLOCKS that tile it.
//   { "pfile_catalog": 2, "dataset": "...",
//     "axes": [
//       {"name":"variants",
//        "keys":[{"name":"CHROM","type":"enum"},
//                {"name":"CHUNK","type":"ordinal_interval","over":"POS"}],
//        "blocks":[{"keys":{"CHROM":"1","CHUNK":"001"},"start":0,"stop":2,
//                   "bounds":{"CHUNK":{"min":10000,"max":20000}}}, ... ]},
//       {"name":"samples",
//        "keys":[{"name":"GROUP","type":"enum"},
//                {"name":"SLICE","type":"ordinal_interval","over":"sample_index"},
//                {"name":"age","type":"value_bin","over":"age"}],
//        "blocks":[{"keys":{"GROUP":"A","SLICE":"0","age":"young"},"start":0,"stop":1,
//                   "bounds":{"age":{"lo":20,"hi":40}}}, ... ]} ],
//     "holes": [{"CHROM":"2","GROUP":"C"}],     // optional: absent cells (see below)
//     "template": {"pfile":"c_{CHROM}_{CHUNK}_{GROUP}_{SLICE}",
//                  "pvar":"v_{CHROM}_{CHUNK}.pvar", "psam":"s_{GROUP}_{SLICE}.psam"} }
//
// Boundary law (§6): a block's ordinal span is HALF-OPEN `[start, stop)`; the blocks of an axis must
// tile `[0, length)` with no gap or overlap, listed in ordinal order. An `ordinal_interval` over POS
// carries INCLUSIVE `{min,max}` bounds per block; a `value_bin` carries HALF-OPEN `{lo,hi}` bin edges.
// An `ordinal_interval` over the axis ordinal (`sample_index`) needs no bound — its bound IS `[start,stop)`.
//
// Sparse grids: an optional top-level `holes` array lists absent cells. Each entry is an object of
// key->value (string) pairs over declared keys of either axis; it removes every product cell whose
// block labels match ALL its pairs (so `{"CHROM":"2"}` drops the whole CHROM=2 row of the grid). A hole
// that matches no cell is a parse error (a typo must not silently leave the cell present).
//
// Back-compat: the phase-1 single-`enum`-key form (`keys:[{name,type:"enum",values:[{value,start,stop}]}]`,
// no `blocks`) is still accepted and NORMALIZED into one block per value at parse time, so the phase-1
// fixtures stay green.
//
// `{KEY}` in a template is replaced by the key's value; `pfile` must reference at least one key of
// EACH axis (distinct cells would otherwise collide); `pvar` may reference only variants-axis keys,
// `psam` only samples-axis keys (§5 invariant, generalized).
// ---------------------------------------------------------------------------

enum class AxisKeyType { ENUM, ORDINAL_INTERVAL, VALUE_BIN };

//! A key's schema on an axis. `over` names the quantity an interval key ranges over: "POS" or
//! "sample_index" for `ordinal_interval`; the value name (e.g. "age") for `value_bin`; "" for `enum`.
struct AxisKey {
	string name;
	AxisKeyType type = AxisKeyType::ENUM;
	string over;
	//! True when each block must carry explicit [lo,hi] bounds for this key (a `value_bin`, or an
	//! `ordinal_interval` over POS). An `ordinal_interval` over the axis ordinal (sample_index) uses
	//! the block's own `[start,stop)` and declares no separate bound.
	bool needs_bounds = false;
};

//! Typed per-block bounds for one interval key. Semantics follow the key type (set at prune time):
//! an `ordinal_interval` over POS is INCLUSIVE `[lo,hi]`; a `value_bin` is HALF-OPEN `[lo,hi)`.
struct BlockBound {
	string key_name;
	int64_t lo = 0;
	int64_t hi = 0;
};

//! One block tiling an axis: the key->value label for every declared key, its half-open ordinal span
//! `[start,stop)` on the axis, and per-key typed bounds.
struct AxisBlock {
	vector<std::pair<string, string>> key_values; //! ordered, one per declared key
	int64_t ordinal_start = 0;
	int64_t ordinal_stop = 0;
	vector<BlockBound> bounds;

	const BlockBound *FindBound(const string &key_name) const {
		for (auto &b : bounds) {
			if (b.key_name == key_name) {
				return &b;
			}
		}
		return nullptr;
	}
	//! The block's label for a key, or "" when the key is not on this axis.
	string ValueOf(const string &key_name) const {
		for (auto &kv : key_values) {
			if (kv.first == key_name) {
				return kv.second;
			}
		}
		return "";
	}
};

//! An axis (`variants` | `samples`), its ordered keys, and the blocks that tile it.
struct AxisDef {
	string name;
	vector<AxisKey> keys;     //! ordered outer->inner
	vector<AxisBlock> blocks; //! ordered, tile [0, length)

	const AxisKey *FindKey(const string &name) const {
		for (auto &k : keys) {
			if (k.name == name) {
				return &k;
			}
		}
		return nullptr;
	}
	//! The genomic-contig key: the enum key marked `"over":"CHROM"`, else an enum key literally named
	//! `CHROM`, else nullptr. (Parse rejects two marked keys.)
	const AxisKey *ContigKey() const {
		for (auto &k : keys) {
			if (k.type == AxisKeyType::ENUM && k.over == "CHROM") {
				return &k;
			}
		}
		auto k = FindKey("CHROM");
		return (k && k->type == AxisKeyType::ENUM) ? k : nullptr;
	}
};

//! One grid cell, resolved. `keys` is (key name, value) for every variants key then every samples key.
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
	int64_t variant_block = 0; //! index into MatrixCatalog::variants.blocks (for bound lookups)
	int64_t sample_block = 0;  //! index into MatrixCatalog::samples.blocks
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
	//! Absent cells (§3.2 sparse grids). Each hole is a partial key->value map over either axis's keys; it
	//! removes every cell of the variant-block x sample-block product whose labels match ALL its entries.
	vector<vector<std::pair<string, string>>> holes;
};

//! True when `content` is a JSON object whose integer `pfile_catalog` is 2 (the version router).
bool IsMatrixCatalog(const string &content);

//! Parse + validate a v2 matrix catalog. Throws "read_pfile_catalog: catalog parse error: ...".
MatrixCatalog ParseMatrixCatalog(FileSystem &fs, const string &path, const string &content);

//! Cartesian product of the two axes' blocks (variants-major) MINUS the declared holes; resolves each
//! PRESENT cell's `.pgen` from the template and its shared `.pvar`/`.psam`. A hole is skipped before any
//! path is resolved, so it is never touched. Opens no file. Throws when two cells collide on a path.
vector<MatrixCell> BuildCells(const MatrixCatalog &cat, FileSystem &fs);

//! Per-key requested value sets (the `keys :=` predicate on `enum` keys): (key name, allowed values).
using KeyPredicates = vector<std::pair<string, vector<string>>>;

//! One interval-overlap predicate carried by `keys` for an `ordinal_interval`/`value_bin` key:
//! an INCLUSIVE `[lo,hi]` matched against the block's typed bound (overlap semantics per key type).
struct BoundPredicate {
	string key_name;
	int64_t lo = 0;
	int64_t hi = 0;
};

//! The resolved, typed prune predicates for a matrix catalog. The reader parses the named parameters
//! (region/variant_index/samples/keys) into this struct; PruneCells applies them as pure geometry.
struct MatrixPredicates {
	KeyPredicates enum_keys;           //! `keys` entries on `enum` keys — membership
	vector<BoundPredicate> bound_keys; //! `keys` entries on `ordinal_interval`/`value_bin` keys — overlap

	bool has_region = false; //! `region` — variant-axis CHROM membership + POS-window overlap
	string region_contig;    //! the contig matched VERBATIM against the variant-axis `CHROM` enum key
	int64_t region_pos_lo = std::numeric_limits<int64_t>::lowest(); //! inclusive POS window (bare contig -> all)
	int64_t region_pos_hi = std::numeric_limits<int64_t>::max();

	bool has_variant_index = false; //! `variant_index := [lo,hi]` — variant-axis ordinal span overlap (inclusive)
	int64_t vi_lo = 0;
	int64_t vi_hi = 0;

	bool has_samples = false; //! `samples := LIST(INTEGER)` — sample-axis ordinal point membership
	vector<int64_t> sample_points;
};

//! Prune: keep a cell iff it passes EVERY predicate (product across axes, §6). `enum`/bound `keys`
//! prune whichever axis the key lives on; `region`/`variant_index` prune the variant axis; `samples`
//! prunes the sample axis. An empty result is not an error.
vector<MatrixCell> PruneCells(const MatrixCatalog &cat, const vector<MatrixCell> &cells, const MatrixPredicates &preds);

} // namespace duckdb
