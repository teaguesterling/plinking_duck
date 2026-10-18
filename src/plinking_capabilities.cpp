#include "plinking_capabilities.hpp"
#include "register_helper.hpp"
#include "duckdb_compat.hpp" // DUCKDB_SCALAR_BIND_*, CompatExprBindInfo, CompatReferenceConstant, CompatMakeName
#include "named_parameter_compat.hpp" // TableFunctionDeclaresParam, CompatDerefTableFunction

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
// exact reported set and proves each reported token binds. Note the asymmetry:
// the test catches a reader LOSING a curated param (the token drops, set-
// equality fails) — but by construction it cannot detect a reader GAINING a
// param that is not in the curated lists. That is an accepted under-claim: an
// uncurated param is simply not a capability token, so there is nothing to
// catch. Surfacing a newly-meaningful param is a deliberate edit here, not an
// automatic one.

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
	// v2.0 dropped the typed Catalog::GetEntry<T>(catalog, schema, name, ...) template.
	// Use the non-throwing instance form on the system catalog (where extension table
	// functions live) and check the type ourselves — the type-arg form does not filter
	// by type reliably (see duckdb_func_apply's note). CompatMakeName gives a string on
	// v1.5 / an Identifier on v2.0 for the name argument.
	auto &system_catalog = Catalog::GetSystemCatalog(context);
	auto entry = system_catalog.GetEntry(context, CatalogType::TABLE_FUNCTION_ENTRY, DEFAULT_SCHEMA,
	                                     CompatMakeName(name), OnEntryNotFound::RETURN_NULL);
	if (!entry || entry->type != CatalogType::TABLE_FUNCTION_ENTRY) {
		return nullptr;
	}
	return &entry->Cast<TableFunctionCatalogEntry>();
}

// True iff `reader` is registered AND EVERY overload exposes a named parameter
// named `param`. A capability claim means "whichever overload the consumer
// calls, this param works", so the param must be present on all overloads (not
// just one) — an asymmetric registration under-claims rather than over-claims.
// A missing reader (no overloads) drops the token (returns false), never throws.
bool ReaderHasParam(ClientContext &context, const string &reader, const string &param) {
	auto entry = LookupTableFunction(context, reader);
	if (!entry || entry->functions.functions.empty()) {
		return false;
	}
	for (auto &fn : entry->functions.functions) {
		// v2.0 holds declared params as typed kwargs on the signature (not a flat map),
		// and the overload list is shared_ptr<const TableFunction> rather than by value.
		if (!TableFunctionDeclaresParam(CompatDerefTableFunction(fn), param)) {
			return false;
		}
	}
	return true;
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

// v2.0 collapsed the scalar bind signature into BindScalarFunctionInput& — see the
// DUCKDB_SCALAR_BIND_* macros in duckdb_compat.hpp.
unique_ptr<FunctionData> CapabilitiesBind(DUCKDB_SCALAR_BIND_PARAMS) {
	return make_uniq<CapabilitiesBindData>(ComputeCapabilities(DUCKDB_SCALAR_BIND_CONTEXT));
}

void CapabilitiesFunc(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
	// v2.0 made BoundFunctionExpression::bind_info private (reached via BindInfo()).
	auto &bind_data = CompatExprBindInfo(func_expr).Cast<CapabilitiesBindData>();

	vector<Value> elems;
	elems.reserve(bind_data.tokens.size());
	for (auto &tok : bind_data.tokens) {
		elems.emplace_back(tok);
	}
	// The result is constant for the whole build; broadcast one constant value.
	// v2.0's Vector::Reference(Value) requires an explicit row count.
	CompatReferenceConstant(result, Value::LIST(LogicalType::VARCHAR, std::move(elems)), args.size());
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
