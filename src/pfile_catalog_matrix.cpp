#include "pfile_catalog_matrix.hpp"
#include "pfile_catalog_util.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"

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

AxisDef ParseAxis(const string &path, yyjson_val *axis, const char *expected_name) {
	AxisDef out;
	out.name = expected_name;
	yyjson_val *keys = yyjson_obj_get(axis, "keys");
	if (!yyjson_is_arr(keys) || yyjson_arr_size(keys) != 1) {
		ParseError(path, StringUtil::Format("axis '%s' needs a 'keys' array with exactly one key (this reader "
		                                    "supports one key per axis)",
		                                    expected_name));
	}
	yyjson_val *k = yyjson_arr_get_first(keys);
	yyjson_val *kname = yyjson_is_obj(k) ? yyjson_obj_get(k, "name") : nullptr;
	if (!kname || !yyjson_is_str(kname)) {
		ParseError(path, StringUtil::Format("axis '%s' key needs a string 'name'", expected_name));
	}
	out.key.name = yyjson_get_str(kname);
	yyjson_val *ktype = yyjson_obj_get(k, "type");
	if (!ktype || !yyjson_is_str(ktype)) {
		ParseError(path, StringUtil::Format("key '%s' needs a string 'type'", out.key.name));
	}
	out.key.type = yyjson_get_str(ktype);
	if (out.key.type != "enum") {
		ParseError(path, StringUtil::Format("key '%s' has unknown key type '%s' (supported: 'enum')", out.key.name,
		                                    out.key.type));
	}
	yyjson_val *values = yyjson_obj_get(k, "values");
	if (!yyjson_is_arr(values) || yyjson_arr_size(values) == 0) {
		ParseError(path, StringUtil::Format("key '%s' needs a non-empty 'values' array", out.key.name));
	}
	size_t idx, max;
	yyjson_val *v;
	int64_t expect_start = 0;
	yyjson_arr_foreach(values, idx, max, v) {
		string where = StringUtil::Format("key '%s' values[%llu]", out.key.name, (unsigned long long)idx);
		yyjson_val *val = yyjson_is_obj(v) ? yyjson_obj_get(v, "value") : nullptr;
		AxisValue av;
		if (!val || !yyjson_is_str(val)) {
			ParseError(path, where + " needs a string 'value'");
		}
		av.value = yyjson_get_str(val);
		if (!JsonGetInt64(yyjson_obj_get(v, "start"), av.start) || !JsonGetInt64(yyjson_obj_get(v, "stop"), av.stop)) {
			ParseError(path, where + " needs integer 'start' and 'stop'");
		}
		if (av.stop <= av.start) {
			ParseError(path, StringUtil::Format("%s span [%lld, %lld) is empty or inverted (stop must be > start, §6)",
			                                    where, (long long)av.start, (long long)av.stop));
		}
		if (av.start != expect_start) {
			ParseError(path, StringUtil::Format("%s span starts at %lld but the axis tiling expects %lld (blocks "
			                                    "must tile the axis with no gap or overlap, §6)",
			                                    where, (long long)av.start, (long long)expect_start));
		}
		expect_start = av.stop;
		for (auto &prev : out.key.values) {
			if (prev.value == av.value) {
				ParseError(path, StringUtil::Format("key '%s' declares value '%s' twice", out.key.name, av.value));
			}
		}
		out.key.values.push_back(std::move(av));
	}
	return out;
}

//! Every `{NAME}` in `tmpl` must be a key in `allowed`; a key of the other axis is the wrong-axis
//! error (§5 invariant, generalized), an unknown name is an undeclared-key error.
void ValidateTemplate(const string &path, const char *which, const string &tmpl, const AxisKey *allowed,
                      const AxisKey *other) {
	size_t i = 0;
	bool uses_own = false;
	while ((i = tmpl.find('{', i)) != string::npos) {
		auto close = tmpl.find('}', i);
		if (close == string::npos) {
			ParseError(path, StringUtil::Format("template.%s '%s' has an unclosed '{'", which, tmpl));
		}
		string name = tmpl.substr(i + 1, close - i - 1);
		if (allowed && name == allowed->name) {
			uses_own = true;
		} else if (other && name == other->name) {
			ParseError(path, StringUtil::Format("template.%s '%s' references key '%s' of the other axis — a "
			                                    "companion is shared only along its own axis (spec §5)",
			                                    which, tmpl, name));
		} else {
			ParseError(path, StringUtil::Format("template.%s '%s' references undeclared key '%s'", which, tmpl, name));
		}
		i = close + 1;
	}
	if (!uses_own) {
		ParseError(path, StringUtil::Format("template.%s '%s' does not reference its own axis's key '%s' — it would "
		                                    "share one companion across every block of that axis (spec §3.1/§5)",
		                                    which, tmpl, allowed->name));
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
	if (cat.variants.key.name == cat.samples.key.name) {
		ParseError(path, StringUtil::Format("both axes use the key name '%s'", cat.variants.key.name));
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

	// The .pgen path may use either axis's key; each companion only its own axis's (§3.1 / §5).
	bool uses_variants = false, uses_samples = false;
	for (size_t pos = 0; (pos = cat.pfile_template.find('{', pos)) != string::npos;) {
		auto close = cat.pfile_template.find('}', pos);
		string name = close == string::npos ? "" : cat.pfile_template.substr(pos + 1, close - pos - 1);
		if (close == string::npos || (name != cat.variants.key.name && name != cat.samples.key.name)) {
			ParseError(path, StringUtil::Format("template.pfile '%s' references undeclared key '%s'",
			                                    cat.pfile_template, name));
		}
		uses_variants |= name == cat.variants.key.name;
		uses_samples |= name == cat.samples.key.name;
		pos = close + 1;
	}
	if (!uses_variants || !uses_samples) {
		ParseError(path, StringUtil::Format("template.pfile '%s' must reference the key of every axis ('%s' and "
		                                    "'%s'), or distinct cells would resolve to the same .pgen",
		                                    cat.pfile_template, cat.variants.key.name, cat.samples.key.name));
	}
	ValidateTemplate(path, "pvar", cat.pvar_template, &cat.variants.key, &cat.samples.key);
	ValidateTemplate(path, "psam", cat.psam_template, &cat.samples.key, &cat.variants.key);
	return cat;
}

vector<MatrixCell> BuildCells(const MatrixCatalog &cat, FileSystem &fs) {
	vector<MatrixCell> cells;
	int64_t ordinal = 0;
	for (auto &v : cat.variants.key.values) {
		for (auto &s : cat.samples.key.values) {
			MatrixCell c;
			c.keys = {{cat.variants.key.name, v.value}, {cat.samples.key.name, s.value}};
			c.prefix = ResolvePath(fs, cat.base_dir, Substitute(cat.pfile_template, c.keys), true);
			c.pvar = ResolvePath(fs, cat.base_dir, Substitute(cat.pvar_template, c.keys), true);
			c.psam = ResolvePath(fs, cat.base_dir, Substitute(cat.psam_template, c.keys), true);
			c.variant_start = v.start;
			c.variant_stop = v.stop;
			c.sample_start = s.start;
			c.sample_stop = s.stop;
			c.source_ordinal = ordinal++;
			cells.push_back(std::move(c));
		}
	}
	return cells;
}

} // namespace duckdb
