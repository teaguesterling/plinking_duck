#include "pfile_catalog_matrix.hpp"
#include "pfile_catalog_util.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"

#include <unordered_set>

using namespace duckdb_yyjson; // NOLINT

namespace duckdb {

// ===========================================================================
// v2 matrix catalog: parse + dense cell manifest. Error prefixes mirror the 1-D reader
// ("read_pfile_catalog: catalog parse error: ...", InvalidInputException).
// ===========================================================================

namespace {

[[noreturn]] void ParseError(const string &path, const string &what) {
	throw InvalidInputException("read_pfile_catalog: catalog parse error: '%s' %s", path, what);
}

//! Parse a key's schema entry `{name, type, over?}`. `over` is required for interval keys and must
//! name the quantity the catalog ranges that key over ("POS" / "sample_index" for an axis ordinal).
AxisKey ParseKeySchema(const string &path, yyjson_val *k, const char *axis_name) {
	AxisKey key;
	yyjson_val *kname = yyjson_is_obj(k) ? yyjson_obj_get(k, "name") : nullptr;
	if (!kname || !yyjson_is_str(kname)) {
		ParseError(path, StringUtil::Format("axis '%s' key needs a string 'name'", axis_name));
	}
	key.name = yyjson_get_str(kname);
	yyjson_val *ktype = yyjson_obj_get(k, "type");
	if (!ktype || !yyjson_is_str(ktype)) {
		ParseError(path, StringUtil::Format("key '%s' needs a string 'type'", key.name));
	}
	string type = yyjson_get_str(ktype);
	yyjson_val *over = yyjson_obj_get(k, "over");
	bool is_variants = string(axis_name) == "variants";
	if (type == "enum") {
		key.type = AxisKeyType::ENUM;
		key.needs_bounds = false;
		// An enum key may mark itself the genomic-contig dimension with "over":"CHROM" (variants axis only).
		if (over) {
			if (!yyjson_is_str(over) || string(yyjson_get_str(over)) != "CHROM") {
				ParseError(path, StringUtil::Format("key '%s' (enum) may only declare 'over':\"CHROM\"", key.name));
			}
			if (!is_variants) {
				ParseError(path, StringUtil::Format("key '%s' on the '%s' axis cannot be 'over':\"CHROM\" (the "
				                                    "genomic contig is a variants-axis dimension)",
				                                    key.name, axis_name));
			}
			key.over = "CHROM";
		}
	} else if (type == "ordinal_interval") {
		key.type = AxisKeyType::ORDINAL_INTERVAL;
		if (!over || !yyjson_is_str(over)) {
			ParseError(path, StringUtil::Format("key '%s' (ordinal_interval) needs a string 'over' "
			                                    "(the axis's ordinal 'POS'/'sample_index')",
			                                    key.name));
		}
		key.over = yyjson_get_str(over);
		const char *axis_ordinal = is_variants ? "POS" : "sample_index";
		if (key.over != axis_ordinal) {
			ParseError(path, StringUtil::Format("key '%s' (ordinal_interval) is on the '%s' axis and must range "
			                                    "'over' '%s', not '%s'",
			                                    key.name, axis_name, axis_ordinal, key.over));
		}
		// An ordinal_interval over POS needs explicit per-block POS bounds; over the axis ordinal, the
		// block's own [start,stop) IS the bound.
		key.needs_bounds = (key.over == "POS");
	} else if (type == "value_bin") {
		key.type = AxisKeyType::VALUE_BIN;
		if (!over || !yyjson_is_str(over)) {
			ParseError(path,
			           StringUtil::Format("key '%s' (value_bin) needs a string 'over' (the binned value)", key.name));
		}
		key.over = yyjson_get_str(over);
		key.needs_bounds = true;
	} else {
		ParseError(path, StringUtil::Format("key '%s' has unknown key type '%s' (supported: 'enum', "
		                                    "'ordinal_interval', 'value_bin')",
		                                    key.name, type));
	}
	return key;
}

//! The phase-1 single-`enum`-key form: one block per `values[{value,start,stop}]` entry. Kept so the
//! phase-1 fixtures parse unchanged. The §6 tiling law (contiguous half-open spans, no dup value) is
//! enforced here exactly as before.
void NormalizeLegacyValues(const string &path, yyjson_val *k, AxisDef &axis, const char *axis_name) {
	yyjson_val *values = yyjson_obj_get(k, "values");
	if (!yyjson_is_arr(values) || yyjson_arr_size(values) == 0) {
		ParseError(path, StringUtil::Format("key '%s' needs a non-empty 'values' array", axis.keys[0].name));
	}
	if (axis.keys[0].type != AxisKeyType::ENUM) {
		ParseError(path, StringUtil::Format("key '%s' uses the single-key 'values' form, which is only for an "
		                                    "'enum' key; declare 'blocks' instead",
		                                    axis.keys[0].name));
	}
	const string &kname = axis.keys[0].name;
	size_t idx, max;
	yyjson_val *v;
	int64_t expect_start = 0;
	yyjson_arr_foreach(values, idx, max, v) {
		string where = StringUtil::Format("key '%s' values[%llu]", kname, (unsigned long long)idx);
		yyjson_val *val = yyjson_is_obj(v) ? yyjson_obj_get(v, "value") : nullptr;
		if (!val || !yyjson_is_str(val)) {
			ParseError(path, where + " needs a string 'value'");
		}
		AxisBlock block;
		string value = yyjson_get_str(val);
		if (!JsonGetInt64(yyjson_obj_get(v, "start"), block.ordinal_start) ||
		    !JsonGetInt64(yyjson_obj_get(v, "stop"), block.ordinal_stop)) {
			ParseError(path, where + " needs integer 'start' and 'stop'");
		}
		if (block.ordinal_stop <= block.ordinal_start) {
			ParseError(path, StringUtil::Format("%s span [%lld, %lld) is empty or inverted (stop must be > start, §6)",
			                                    where, (long long)block.ordinal_start, (long long)block.ordinal_stop));
		}
		if (block.ordinal_start != expect_start) {
			ParseError(path, StringUtil::Format("%s span starts at %lld but the axis tiling expects %lld (blocks "
			                                    "must tile the axis with no gap or overlap, §6)",
			                                    where, (long long)block.ordinal_start, (long long)expect_start));
		}
		expect_start = block.ordinal_stop;
		for (auto &prev : axis.blocks) {
			if (prev.key_values[0].second == value) {
				ParseError(path, StringUtil::Format("key '%s' declares value '%s' twice", kname, value));
			}
		}
		block.key_values = {{kname, value}};
		axis.blocks.push_back(std::move(block));
	}
	(void)axis_name;
}

//! The `blocks[]` tiling form: each block lists a value for EVERY declared key, its half-open ordinal
//! span, and per-key typed bounds (required for a `value_bin` or an `ordinal_interval` over POS).
void ParseBlocks(const string &path, yyjson_val *blocks, AxisDef &axis, const char *axis_name) {
	if (!yyjson_is_arr(blocks) || yyjson_arr_size(blocks) == 0) {
		ParseError(path, StringUtil::Format("axis '%s' needs a non-empty 'blocks' array", axis_name));
	}
	size_t idx, max;
	yyjson_val *b;
	int64_t expect_start = 0;
	std::unordered_set<string> seen_combos;
	yyjson_arr_foreach(blocks, idx, max, b) {
		string where = StringUtil::Format("axis '%s' blocks[%llu]", axis_name, (unsigned long long)idx);
		if (!yyjson_is_obj(b)) {
			ParseError(path, where + " must be an object {keys, start, stop[, bounds]}");
		}
		AxisBlock block;
		if (!JsonGetInt64(yyjson_obj_get(b, "start"), block.ordinal_start) ||
		    !JsonGetInt64(yyjson_obj_get(b, "stop"), block.ordinal_stop)) {
			ParseError(path, where + " needs integer 'start' and 'stop'");
		}
		if (block.ordinal_stop <= block.ordinal_start) {
			ParseError(path, StringUtil::Format("%s span [%lld, %lld) is empty or inverted (stop must be > start, §6)",
			                                    where, (long long)block.ordinal_start, (long long)block.ordinal_stop));
		}
		if (block.ordinal_start != expect_start) {
			ParseError(path, StringUtil::Format("%s span starts at %lld but the axis tiling expects %lld (blocks "
			                                    "must tile the axis with no gap or overlap, §6)",
			                                    where, (long long)block.ordinal_start, (long long)expect_start));
		}
		expect_start = block.ordinal_stop;

		yyjson_val *kv = yyjson_obj_get(b, "keys");
		if (!yyjson_is_obj(kv)) {
			ParseError(path, where + " needs a 'keys' object {<key>:<value>, ...}");
		}
		yyjson_val *bounds = yyjson_obj_get(b, "bounds");
		if (bounds && !yyjson_is_obj(bounds)) {
			ParseError(path, where + " 'bounds' must be an object");
		}
		string combo;
		for (auto &key : axis.keys) {
			yyjson_val *kval = yyjson_obj_get(kv, key.name.c_str());
			if (!kval || !yyjson_is_str(kval)) {
				ParseError(path, StringUtil::Format("%s keys is missing a string value for declared key '%s'", where,
				                                    key.name));
			}
			string value = yyjson_get_str(kval);
			block.key_values.emplace_back(key.name, value);
			combo += value;
			combo += '\x01';

			yyjson_val *kb = bounds ? yyjson_obj_get(bounds, key.name.c_str()) : nullptr;
			if (key.needs_bounds) {
				if (!yyjson_is_obj(kb)) {
					const char *lohi = key.type == AxisKeyType::VALUE_BIN ? "{lo,hi}" : "{min,max}";
					ParseError(path,
					           StringUtil::Format("%s needs bounds.%s %s for the %s key '%s'", where, key.name, lohi,
					                              key.type == AxisKeyType::VALUE_BIN ? "value_bin" : "ordinal_interval",
					                              key.name));
				}
				BlockBound bound;
				bound.key_name = key.name;
				// value_bin carries half-open {lo,hi}; ordinal_interval over POS carries inclusive {min,max}.
				const char *lo_field = key.type == AxisKeyType::VALUE_BIN ? "lo" : "min";
				const char *hi_field = key.type == AxisKeyType::VALUE_BIN ? "hi" : "max";
				if (!JsonGetInt64(yyjson_obj_get(kb, lo_field), bound.lo) ||
				    !JsonGetInt64(yyjson_obj_get(kb, hi_field), bound.hi)) {
					ParseError(path, StringUtil::Format("%s bounds.%s needs integer '%s' and '%s'", where, key.name,
					                                    lo_field, hi_field));
				}
				if (bound.hi < bound.lo) {
					ParseError(path, StringUtil::Format("%s bounds.%s [%lld, %lld] is inverted (hi/max must be >= "
					                                    "lo/min)",
					                                    where, key.name, (long long)bound.lo, (long long)bound.hi));
				}
				block.bounds.push_back(std::move(bound));
			}
		}
		if (!seen_combos.insert(combo).second) {
			ParseError(path, StringUtil::Format("%s declares a key-value combination twice", where));
		}
		axis.blocks.push_back(std::move(block));
	}
}

AxisDef ParseAxis(const string &path, yyjson_val *axis_val, const char *expected_name) {
	AxisDef out;
	out.name = expected_name;
	yyjson_val *keys = yyjson_obj_get(axis_val, "keys");
	if (!yyjson_is_arr(keys) || yyjson_arr_size(keys) == 0) {
		ParseError(path, StringUtil::Format("axis '%s' needs a non-empty 'keys' array", expected_name));
	}
	size_t idx, max;
	yyjson_val *k;
	yyjson_arr_foreach(keys, idx, max, k) {
		AxisKey key = ParseKeySchema(path, k, expected_name);
		for (auto &prev : out.keys) {
			if (prev.name == key.name) {
				ParseError(path, StringUtil::Format("axis '%s' declares key '%s' twice", expected_name, key.name));
			}
		}
		if (key.type == AxisKeyType::ENUM && key.over == "CHROM") {
			for (auto &prev : out.keys) {
				if (prev.type == AxisKeyType::ENUM && prev.over == "CHROM") {
					ParseError(path, StringUtil::Format("axis '%s' keys '%s' and '%s' both claim 'over':\"CHROM\"; "
					                                    "at most one key may be the contig dimension",
					                                    expected_name, prev.name, key.name));
				}
			}
		}
		out.keys.push_back(std::move(key));
	}

	yyjson_val *blocks = yyjson_obj_get(axis_val, "blocks");
	bool has_legacy_values = out.keys.size() == 1 && yyjson_obj_get(yyjson_arr_get_first(keys), "values");
	if (out.keys.size() > 1 && blocks) {
		// 'values' is the single-key form only; beside 'blocks' on a multi-key axis it is stale/ambiguous, never
		// silently ignored. (Without 'blocks' the generic needs-a-'blocks' error below applies.)
		size_t kidx, kmax;
		yyjson_val *kk;
		yyjson_arr_foreach(keys, kidx, kmax, kk) {
			if (yyjson_obj_get(kk, "values")) {
				ParseError(path, StringUtil::Format("axis '%s' key '%s' carries a 'values' array, which is only valid "
				                                    "on a single-key axis; declare the tiling in 'blocks'",
				                                    expected_name, out.keys[kidx].name));
			}
		}
	}
	if (blocks) {
		if (has_legacy_values) {
			ParseError(path, StringUtil::Format("axis '%s' has both 'blocks' and a single-key 'values' list; use "
			                                    "one or the other",
			                                    expected_name));
		}
		ParseBlocks(path, blocks, out, expected_name);
	} else if (out.keys.size() == 1 && has_legacy_values) {
		NormalizeLegacyValues(path, yyjson_arr_get_first(keys), out, expected_name);
	} else {
		ParseError(path, StringUtil::Format("axis '%s' with more than one key, or without a single-key 'values' "
		                                    "list, needs a 'blocks' array (the tiling)",
		                                    expected_name));
	}
	return out;
}

//! Every `{NAME}` in `tmpl` must be a key of `allowed`; a key of `other` is the wrong-axis error
//! (§5 invariant, generalized), an unknown name is an undeclared-key error. Requires at least one
//! `allowed` key (a companion shared across every block of its axis would be wrong, §3.1/§5).
void ValidateTemplate(const string &path, const char *which, const string &tmpl, const AxisDef &allowed,
                      const AxisDef &other) {
	size_t i = 0;
	bool uses_own = false;
	while ((i = tmpl.find('{', i)) != string::npos) {
		auto close = tmpl.find('}', i);
		if (close == string::npos) {
			ParseError(path, StringUtil::Format("template.%s '%s' has an unclosed '{'", which, tmpl));
		}
		string name = tmpl.substr(i + 1, close - i - 1);
		if (allowed.FindKey(name)) {
			uses_own = true;
		} else if (other.FindKey(name)) {
			ParseError(path, StringUtil::Format("template.%s '%s' references key '%s' of the other axis — a "
			                                    "companion is shared only along its own axis (spec §5)",
			                                    which, tmpl, name));
		} else {
			ParseError(path, StringUtil::Format("template.%s '%s' references undeclared key '%s'", which, tmpl, name));
		}
		i = close + 1;
	}
	if (!uses_own) {
		ParseError(path, StringUtil::Format("template.%s '%s' does not reference its own axis's key(s) — it would "
		                                    "share one companion across every block of that axis (spec §3.1/§5)",
		                                    which, tmpl));
	}
}

string Substitute(const string &tmpl, const vector<std::pair<string, string>> &keys) {
	string out = tmpl;
	for (auto &kv : keys) {
		out = StringUtil::Replace(out, "{" + kv.first + "}", kv.second);
	}
	return out;
}

} // namespace

bool IsMatrixCatalog(const string &content) {
	YyDocGuard guard;
	guard.doc = yyjson_read(content.c_str(), content.size(), 0);
	if (!guard.doc) {
		return false;
	}
	yyjson_val *root = yyjson_doc_get_root(guard.doc);
	int64_t version = 0;
	return yyjson_is_obj(root) && JsonGetInt64(yyjson_obj_get(root, "pfile_catalog"), version) && version == 2;
}

namespace {

//! True when the hole's every (key, value) matches the labels of this variant-block x sample-block pair.
bool HoleMatches(const MatrixCatalog &cat, const vector<std::pair<string, string>> &hole, const AxisBlock &vb,
                 const AxisBlock &sb) {
	for (auto &kv : hole) {
		auto &block = cat.variants.FindKey(kv.first) ? vb : sb;
		if (block.ValueOf(kv.first) != kv.second) {
			return false;
		}
	}
	return true;
}

bool IsHole(const MatrixCatalog &cat, const AxisBlock &vb, const AxisBlock &sb) {
	for (auto &hole : cat.holes) {
		if (HoleMatches(cat, hole, vb, sb)) {
			return true;
		}
	}
	return false;
}

void ParseHoles(const string &path, yyjson_val *root, MatrixCatalog &cat) {
	yyjson_val *holes = yyjson_obj_get(root, "holes");
	if (!holes) {
		return;
	}
	if (!yyjson_is_arr(holes)) {
		ParseError(path, "'holes' must be an array of {key: value} objects");
	}
	size_t idx, max;
	yyjson_val *h;
	yyjson_arr_foreach(holes, idx, max, h) {
		if (!yyjson_is_obj(h) || yyjson_obj_size(h) == 0) {
			ParseError(path, StringUtil::Format("holes[%llu] must be a non-empty object of key: value pairs",
			                                    (unsigned long long)idx));
		}
		vector<std::pair<string, string>> hole;
		size_t i2, m2;
		yyjson_val *k, *v;
		yyjson_obj_foreach(h, i2, m2, k, v) {
			string key = yyjson_get_str(k);
			if (!cat.variants.FindKey(key) && !cat.samples.FindKey(key)) {
				ParseError(path,
				           StringUtil::Format("holes[%llu] names undeclared key '%s'", (unsigned long long)idx, key));
			}
			if (!yyjson_is_str(v)) {
				ParseError(path,
				           StringUtil::Format("holes[%llu]['%s'] must be a string", (unsigned long long)idx, key));
			}
			hole.emplace_back(key, yyjson_get_str(v));
		}
		bool matches_any = false;
		for (auto &vb : cat.variants.blocks) {
			for (auto &sb : cat.samples.blocks) {
				matches_any |= HoleMatches(cat, hole, vb, sb);
			}
		}
		if (!matches_any) {
			ParseError(path, StringUtil::Format("holes[%llu] matches no cell of the grid", (unsigned long long)idx));
		}
		cat.holes.push_back(std::move(hole));
	}
}

} // namespace

MatrixCatalog ParseMatrixCatalog(FileSystem &fs, const string &path, const string &content) {
	YyDocGuard guard;
	guard.doc = yyjson_read(content.c_str(), content.size(), 0);
	if (!guard.doc) {
		ParseError(path, "is not valid JSON");
	}
	yyjson_val *root = yyjson_doc_get_root(guard.doc);
	if (!yyjson_is_obj(root)) {
		ParseError(path, "top level must be a JSON object");
	}

	MatrixCatalog cat;
	cat.base_dir = DirName(AbsPath(fs, path));
	if (!JsonGetInt64(yyjson_obj_get(root, "pfile_catalog"), cat.version) || cat.version != 2) {
		ParseError(path, "is not a version-2 matrix catalog");
	}

	yyjson_val *ds = yyjson_obj_get(root, "dataset");
	if (!ds || !yyjson_is_str(ds)) {
		ParseError(path, "missing required string field 'dataset'");
	}
	cat.dataset = yyjson_get_str(ds);

	// v1 (1-D / single-fileset) keys cannot coexist with the v2 'axes' block.
	if (yyjson_obj_get(root, "concat_dim") || yyjson_obj_get(root, "files") || yyjson_obj_get(root, "shared") ||
	    yyjson_obj_get(root, "pfile")) {
		ParseError(path, "has both a v2 'axes' block and v1 ('concat_dim'/'files'/'shared'/'pfile') keys; a "
		                 "catalog is one or the other");
	}

	yyjson_val *axes = yyjson_obj_get(root, "axes");
	if (!yyjson_is_arr(axes)) {
		ParseError(path, "version 2 needs an 'axes' array");
	}
	bool have_variants = false, have_samples = false;
	size_t idx, max;
	yyjson_val *ax;
	yyjson_arr_foreach(axes, idx, max, ax) {
		yyjson_val *name = yyjson_is_obj(ax) ? yyjson_obj_get(ax, "name") : nullptr;
		string n = name && yyjson_is_str(name) ? yyjson_get_str(name) : "";
		if (n == "variants" && !have_variants) {
			cat.variants = ParseAxis(path, ax, "variants");
			have_variants = true;
		} else if (n == "samples" && !have_samples) {
			cat.samples = ParseAxis(path, ax, "samples");
			have_samples = true;
		} else {
			ParseError(path, StringUtil::Format("axes[%llu] must be a unique axis named 'variants' or 'samples'",
			                                    (unsigned long long)idx));
		}
	}
	if (!have_variants || !have_samples) {
		ParseError(path, "'axes' must declare both a 'variants' and a 'samples' axis");
	}
	// No key name may be shared ACROSS the two axes (a template `{NAME}` must resolve to one axis).
	for (auto &vk : cat.variants.keys) {
		if (cat.samples.FindKey(vk.name)) {
			ParseError(path, StringUtil::Format("both axes use the key name '%s'", vk.name));
		}
	}

	yyjson_val *tmpl = yyjson_obj_get(root, "template");
	if (!yyjson_is_obj(tmpl)) {
		ParseError(path, "version 2 needs a 'template' object {pfile, pvar, psam}");
	}
	auto get_tmpl = [&](const char *field) {
		yyjson_val *t = yyjson_obj_get(tmpl, field);
		if (!t || !yyjson_is_str(t)) {
			ParseError(path, StringUtil::Format("template needs a string '%s'", field));
		}
		return string(yyjson_get_str(t));
	};
	cat.pfile_template = get_tmpl("pfile");
	cat.pvar_template = get_tmpl("pvar");
	cat.psam_template = get_tmpl("psam");

	// The .pgen path may use either axis's keys; it must reference at least one key of EACH axis, or
	// distinct cells collapse to the same .pgen (§3.1 / §5).
	bool uses_variants = false, uses_samples = false;
	for (size_t pos = 0; (pos = cat.pfile_template.find('{', pos)) != string::npos;) {
		auto close = cat.pfile_template.find('}', pos);
		string name = close == string::npos ? "" : cat.pfile_template.substr(pos + 1, close - pos - 1);
		if (close == string::npos) {
			ParseError(path, StringUtil::Format("template.pfile '%s' has an unclosed '{'", cat.pfile_template));
		}
		if (!cat.variants.FindKey(name) && !cat.samples.FindKey(name)) {
			ParseError(path, StringUtil::Format("template.pfile '%s' references undeclared key '%s'",
			                                    cat.pfile_template, name));
		}
		uses_variants |= cat.variants.FindKey(name) != nullptr;
		uses_samples |= cat.samples.FindKey(name) != nullptr;
		pos = close + 1;
	}
	if (!uses_variants || !uses_samples) {
		ParseError(path, StringUtil::Format("template.pfile '%s' must reference the key of every axis, or distinct "
		                                    "cells would resolve to the same .pgen",
		                                    cat.pfile_template));
	}
	ValidateTemplate(path, "pvar", cat.pvar_template, cat.variants, cat.samples);
	ValidateTemplate(path, "psam", cat.psam_template, cat.samples, cat.variants);
	ParseHoles(path, root, cat);
	return cat;
}

vector<MatrixCell> BuildCells(const MatrixCatalog &cat, FileSystem &fs) {
	vector<MatrixCell> cells;
	std::unordered_set<string> seen_prefix, seen_pvar, seen_psam;
	int64_t ordinal = 0;
	// A companion is registered (and collision-checked) once per block, by the first PRESENT cell of that
	// block — not by grid position, since holes may remove the first cell of a row/column.
	vector<bool> pvar_done(cat.variants.blocks.size(), false), psam_done(cat.samples.blocks.size(), false);
	for (size_t vi = 0; vi < cat.variants.blocks.size(); vi++) {
		auto &vb = cat.variants.blocks[vi];
		for (size_t si = 0; si < cat.samples.blocks.size(); si++) {
			auto &sb = cat.samples.blocks[si];
			if (IsHole(cat, vb, sb)) {
				continue; // absent: never resolve a path for it
			}
			MatrixCell c;
			c.keys = vb.key_values;
			c.keys.insert(c.keys.end(), sb.key_values.begin(), sb.key_values.end());
			c.prefix = ResolvePath(fs, cat.base_dir, Substitute(cat.pfile_template, c.keys), true);
			c.pvar = ResolvePath(fs, cat.base_dir, Substitute(cat.pvar_template, c.keys), true);
			c.psam = ResolvePath(fs, cat.base_dir, Substitute(cat.psam_template, c.keys), true);
			// Each variant block owns one `.pvar`, each sample block one `.psam` (§3.1: an M×N grid has
			// M `.pvar` + N `.psam`). A companion template that omits a block-distinguishing key would pair
			// a block with another block's companion — caught here, not left to coincidental count match.
			if (!pvar_done[vi] && !seen_pvar.insert(c.pvar).second) {
				throw InvalidInputException("read_pfile_catalog: catalog parse error: two variant-blocks resolve to "
				                            "the same .pvar '%s' — the template.pvar omits a variant key that "
				                            "distinguishes them",
				                            c.pvar);
			}
			if (!psam_done[si] && !seen_psam.insert(c.psam).second) {
				throw InvalidInputException("read_pfile_catalog: catalog parse error: two sample-blocks resolve to "
				                            "the same .psam '%s' — the template.psam omits a sample key that "
				                            "distinguishes them",
				                            c.psam);
			}
			if (!seen_prefix.insert(c.prefix).second) {
				throw InvalidInputException("read_pfile_catalog: catalog parse error: two cells resolve to the same "
				                            ".pgen prefix '%s' — the template omits a key that distinguishes them",
				                            c.prefix);
			}
			pvar_done[vi] = psam_done[si] = true;
			c.variant_start = vb.ordinal_start;
			c.variant_stop = vb.ordinal_stop;
			c.sample_start = sb.ordinal_start;
			c.sample_stop = sb.ordinal_stop;
			c.variant_block = (int64_t)vi;
			c.sample_block = (int64_t)si;
			c.source_ordinal = ordinal++;
			cells.push_back(std::move(c));
		}
	}
	return cells;
}

namespace {

//! The interval a bound-key predicate overlaps against, resolved from the cell's block. Returns false
//! when the key is not an interval key on either axis (handled as a predicate error upstream).
bool ResolveBoundInterval(const MatrixCatalog &cat, const MatrixCell &c, const string &key_name, int64_t &lo,
                          int64_t &hi, bool &inclusive) {
	for (auto *axis : {&cat.variants, &cat.samples}) {
		const AxisKey *key = axis->FindKey(key_name);
		if (!key) {
			continue;
		}
		const AxisBlock &block =
		    axis->name == "variants" ? cat.variants.blocks[c.variant_block] : cat.samples.blocks[c.sample_block];
		if (key->type == AxisKeyType::ORDINAL_INTERVAL && key->over != "POS") {
			// over an axis ordinal: the bound IS the block's half-open [start,stop).
			lo = block.ordinal_start;
			hi = block.ordinal_stop;
			inclusive = false;
			return true;
		}
		const BlockBound *b = block.FindBound(key_name);
		if (!b) {
			return false;
		}
		lo = b->lo;
		hi = b->hi;
		inclusive = key->type == AxisKeyType::ORDINAL_INTERVAL; // POS inclusive; value_bin half-open
		return true;
	}
	return false;
}

} // namespace

vector<MatrixCell> PruneCells(const MatrixCatalog &cat, const vector<MatrixCell> &cells,
                              const MatrixPredicates &preds) {
	// Resolve the variant-axis ordinal_interval-over-POS key once (for `region`'s POS window).
	const AxisKey *pos_key = nullptr;
	for (auto &k : cat.variants.keys) {
		if (k.type == AxisKeyType::ORDINAL_INTERVAL && k.over == "POS") {
			pos_key = &k;
			break;
		}
	}

	const AxisKey *contig_key = cat.variants.ContigKey();

	vector<MatrixCell> out;
	for (auto &c : cells) {
		bool keep = true;

		// enum `keys` membership (either axis).
		for (auto &p : preds.enum_keys) {
			for (auto &kv : c.keys) {
				if (kv.first != p.first) {
					continue;
				}
				bool in_set = false;
				for (auto &v : p.second) {
					in_set = in_set || v == kv.second;
				}
				keep = keep && in_set;
			}
		}

		// bound `keys` overlap (ordinal_interval / value_bin, either axis).
		for (auto &bp : preds.bound_keys) {
			int64_t blo, bhi;
			bool inclusive;
			if (!ResolveBoundInterval(cat, c, bp.key_name, blo, bhi, inclusive)) {
				keep = false; // no bound on this block -> cannot overlap
				continue;
			}
			bool overlap = inclusive ? !(bp.hi < blo || bp.lo > bhi)  // inclusive [blo,bhi] vs inclusive [lo,hi]
			                         : (bp.lo < bhi && bp.hi >= blo); // half-open [blo,bhi) vs inclusive [lo,hi]
			keep = keep && overlap;
		}

		// region: variant-axis CHROM membership (verbatim) + POS-window overlap against the POS key.
		if (preds.has_region) {
			string chrom =
			    (contig_key && !c.keys.empty()) ? cat.variants.blocks[c.variant_block].ValueOf(contig_key->name) : "";
			if (chrom != preds.region_contig) {
				keep = false;
			} else if (pos_key) {
				const BlockBound *b = cat.variants.blocks[c.variant_block].FindBound(pos_key->name);
				if (b && (preds.region_pos_hi < b->lo || preds.region_pos_lo > b->hi)) {
					keep = false; // inclusive POS window disjoint from the block's inclusive POS extent
				}
			}
		}

		// variant_index: inclusive [lo,hi] overlaps the variant block's half-open [start,stop).
		if (preds.has_variant_index && !(c.variant_start <= preds.vi_hi && c.variant_stop > preds.vi_lo)) {
			keep = false;
		}

		// samples: any listed index falls in the sample block's half-open [start,stop).
		if (preds.has_samples) {
			bool hit = false;
			for (auto pt : preds.sample_points) {
				if (c.sample_start <= pt && pt < c.sample_stop) {
					hit = true;
					break;
				}
			}
			keep = keep && hit;
		}

		if (keep) {
			out.push_back(c);
		}
	}
	return out;
}

} // namespace duckdb
