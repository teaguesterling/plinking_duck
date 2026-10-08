#include "plinking_capabilities.hpp"
#include "register_helper.hpp"

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_function_catalog_entry.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/enums/catalog_type.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"

namespace duckdb {

// ===========================================================================
// plinking_duck_capabilities() -> VARCHAR[]
// ===========================================================================
//
// Returns the list of capability tokens THIS build of plinking_duck actually
// supports. A consumer detects a feature in pure SQL with, e.g.:
//
//     SELECT list_contains(plinking_duck_capabilities(), 'variant_range');
//
// and detects that this detection mechanism itself is present with:
//
//     SELECT count(*) > 0 FROM duckdb_functions()
//     WHERE function_name = 'plinking_duck_capabilities';
//
// CONTRACT / TOKEN NAMING
// -----------------------
// Tokens are consumer-facing and stable. They are spelled exactly as the
// underlying named-parameter key or table-function name (snake_case), so the
// contract reads identically in SQL and in the registration code.
//
// DERIVED, NOT HAND-MAINTAINED
// ----------------------------
// The reported set is computed from the LIVE catalog at bind time, never from a
// hardcoded literal list of "what we believe this build does". A token is
// emitted ONLY when its underlying registration fact is genuinely present in
// this build. Two kinds of fact:
//
//   * PARAM capability: a named parameter `<token>` is registered on EVERY
//     reader listed for it (we inspect each reader's live named_parameters).
//   * FUNCTION capability: a table function named `<token>` exists in the
//     catalog.
//
// Consequence: renaming/removing a reader param or dropping a function
// automatically drops the token — capabilities() can never claim a feature the
// readers do not actually honor. The curated lists below only decide WHICH
// facts are consumer-meaningful enough to surface as a token (we deliberately
// do not enumerate every internal named parameter); they do not assert that the
// feature works. The companion test (plinking_duck_capabilities.test) pins the
// exact token set in both directions and proves each reported token binds.

namespace {

// A named-parameter capability: `token` is reported iff every reader in
// `readers` exposes a named parameter named `token`.
struct ParamCapability {
	const char *token;
	vector<string> readers;
};

// Curated set of param-backed capabilities (only consumer-meaningful params).
const vector<ParamCapability> &ParamCapabilities() {
	static const vector<ParamCapability> caps = {
	    // variant_range := [lo, hi) — ordinal-range pushdown on read_pfile.
	    {"variant_range", {"read_pfile"}},
	    // with_ordinal := true — appends the __ordinal column; both readers honor it.
	    {"with_ordinal", {"read_pfile", "read_pvar"}},
	};
	return caps;
}

// Curated set of function-name capabilities: reported iff the function exists.
const vector<string> &FunctionCapabilities() {
	static const vector<string> caps = {
	    "pfile_partitions",
	};
	return caps;
}

optional_ptr<TableFunctionCatalogEntry> LookupTableFunction(ClientContext &context, const string &name) {
	return Catalog::GetEntry<TableFunctionCatalogEntry>(context, INVALID_CATALOG, DEFAULT_SCHEMA, name,
	                                                    OnEntryNotFound::RETURN_NULL);
}

// True iff `reader` is registered AND at least one of its overloads exposes a
// named parameter named `param`.
bool ReaderHasParam(ClientContext &context, const string &reader, const string &param) {
	auto entry = LookupTableFunction(context, reader);
	if (!entry) {
		return false;
	}
	for (auto &fn : entry->functions.functions) {
		if (fn.named_parameters.find(param) != fn.named_parameters.end()) {
			return true;
		}
	}
	return false;
}

// Compute the live capability set from the catalog.
vector<string> ComputeCapabilities(ClientContext &context) {
	vector<string> tokens;
	for (auto &cap : ParamCapabilities()) {
		bool present = !cap.readers.empty();
		for (auto &reader : cap.readers) {
			if (!ReaderHasParam(context, reader, cap.token)) {
				present = false;
				break;
			}
		}
		if (present) {
			tokens.emplace_back(cap.token);
		}
	}
	for (auto &fname : FunctionCapabilities()) {
		if (LookupTableFunction(context, fname)) {
			tokens.push_back(fname);
		}
	}
	return tokens;
}

struct CapabilitiesBindData : public FunctionData {
	explicit CapabilitiesBindData(vector<string> tokens_p) : tokens(std::move(tokens_p)) {
	}
	vector<string> tokens;

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<CapabilitiesBindData>(tokens);
	}
	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<CapabilitiesBindData>();
		return tokens == other.tokens;
	}
};

unique_ptr<FunctionData> CapabilitiesBind(ClientContext &context, ScalarFunction &, vector<unique_ptr<Expression>> &) {
	return make_uniq<CapabilitiesBindData>(ComputeCapabilities(context));
}

void CapabilitiesFunc(DataChunk &, ExpressionState &state, Vector &result) {
	auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
	auto &bind_data = func_expr.bind_info->Cast<CapabilitiesBindData>();

	vector<Value> elems;
	elems.reserve(bind_data.tokens.size());
	for (auto &tok : bind_data.tokens) {
		elems.emplace_back(Value(tok));
	}
	// The result is constant for the whole build; broadcast one constant value.
	result.Reference(Value::LIST(LogicalType::VARCHAR, std::move(elems)));
}

} // namespace

void RegisterPlinkingCapabilities(ExtensionLoader &loader) {
	ScalarFunction fn("plinking_duck_capabilities", {}, LogicalType::LIST(LogicalType::VARCHAR), CapabilitiesFunc,
	                  CapabilitiesBind);
	RegisterScalarWithDesc(loader, fn, {},
	                       "List the capability tokens this plinking_duck build supports (for SQL feature detection).",
	                       {"SELECT list_contains(plinking_duck_capabilities(), 'variant_range')"});
}

} // namespace duckdb
