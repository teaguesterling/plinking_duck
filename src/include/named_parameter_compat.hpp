//===----------------------------------------------------------------------===//
//                          plinking_duck
//
// named_parameter_compat.hpp
//
// Declaring and reading a table function's named parameters, on both DuckDB
// lines.
//
//===----------------------------------------------------------------------===//
//
// LINEAGE: lifted verbatim from the fleet-shared named_parameter_compat.hpp
// (origin sitting_duck; reference copy duckdb_func_apply/src/include/). The
// SITTING_DUCK_* macro/group names are kept unchanged so this stays a by-name
// merge with the fleet copies — do not rename them. plinking_duck added nothing
// to this file; if the fleet copy evolves, re-sync rather than fork.
//
// DuckDB v2.0 replaced the flat named-parameter map with a signature carrying
// typed kwargs:
//
//   v1.5   func.named_parameters["peek_size"] = LogicalType::INTEGER;
//   v2.0   func.GetSignature().WithTypedKwargs("options", [](TypedKwargs &o) {
//              o.Add(Identifier("peek_size"), LogicalType::INTEGER);
//          });
//
// and retyped the bind-time map from case_insensitive_map_t<Value> (keyed by
// string) to named_argument_map_t (keyed by Identifier), so `.at("peek_size")`
// no longer compiles: Identifier's constructor from string is explicit by
// design.
//
// WHY #if HERE, WHEN THE REST OF THE EXTENSION USES SFINAE PROBES.
// A flat map and a signature with named kwargs groups are different shapes, and
// the v2.0 spelling names types (TypedKwargs) that do not exist on v1.5 at all.
// A template cannot hide a type that is absent -- non-dependent names are looked
// up when the template is defined, not when it is instantiated. So the branch is
// a preprocessor one, confined to this header.
//
// THE SENTINEL. duckdb/main/capi/capi_function_signature.hpp is absent on v1.5.6
// and arrived on v2.0-cyanoptera in the same work that gave TableFunction its
// signature. Rejected sentinels: identifier.hpp (backported to v1.5, present on
// both) and enums/identifier_case_mode.hpp (v2.0-only but landed BEFORE
// GetSignature at cyanoptera e366461e30 — selected the v2.0 path against a DuckDB
// with no kwargs API, 99 errors). A sentinel must co-vary with the API it gates.
//
#pragma once

#include "duckdb.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/function/table_function.hpp"

#if __has_include("duckdb/main/capi/capi_function_signature.hpp")
#define SITTING_DUCK_HAS_TYPED_KWARGS 1
#endif

namespace duckdb {

//! One named parameter: the name callers write, and the type it accepts.
struct NamedParamSpec {
	const char *name;
	LogicalType type;
};

//! The name of the kwargs group v2.0 collects a table function's options under.
//! v1.5 has no grouping, so this is unused there.
static constexpr const char *SITTING_DUCK_KWARGS_GROUP = "options";

//! The type of TableFunctionBindInput::named_parameters on whichever line we are
//! building against. Spell the map this way in any signature that receives it:
//! naming named_parameter_map_t directly compiles on v1.5 and silently means a
//! DIFFERENT type on v2.0 (which still defines that name, as identifier_map_t,
//! while the bind input actually hands out named_argument_map_t).
#ifdef SITTING_DUCK_HAS_TYPED_KWARGS
using NamedParamMap = named_argument_map_t;
#else
using NamedParamMap = named_parameter_map_t;
#endif

#ifdef SITTING_DUCK_HAS_TYPED_KWARGS

//! Declare a table function's named parameters (v2.0: typed kwargs on the signature).
//! Call this ONCE per function, then Extend for any further parameters: v2.0
//! distinguishes creating the "**kwargs" parameter from adding options to it,
//! and ExtendTypedKwargs throws if the signature has no kwargs parameter yet.
inline void DeclareNamedParameters(TableFunction &func, const vector<NamedParamSpec> &params) {
	func.GetSignature().WithTypedKwargs(Identifier(SITTING_DUCK_KWARGS_GROUP), [&params](TypedKwargs &kwargs) {
		for (const auto &param : params) {
			kwargs.Add(Identifier(string(param.name)), param.type);
		}
	});
}

//! Add more named parameters to a function that has already declared some.
inline void ExtendNamedParameters(TableFunction &func, const vector<NamedParamSpec> &params) {
	func.GetSignature().ExtendTypedKwargs([&params](TypedKwargs &kwargs) {
		for (const auto &param : params) {
			kwargs.Add(Identifier(string(param.name)), param.type);
		}
	});
}

//! The value bound to `name`, or nullptr when the call left it out.
//! v2.0 keys this map by Identifier, and Identifier's string constructor is
//! explicit, so the key has to be built rather than passed as a literal.
inline const Value *FindNamedParameter(const named_argument_map_t &params, const char *name) {
	auto entry = params.find(Identifier(string(name)));
	if (entry == params.end()) {
		return nullptr;
	}
	return &entry->second;
}

#else

//! Declare a table function's named parameters (v1.5: the flat map).
inline void DeclareNamedParameters(TableFunction &func, const vector<NamedParamSpec> &params) {
	for (const auto &param : params) {
		func.named_parameters[param.name] = param.type;
	}
}

//! Add more named parameters to a function that has already declared some. On
//! this line there is no distinction -- the map does not care when a key arrives
//! -- but the call sites keep the same shape on both lines.
inline void ExtendNamedParameters(TableFunction &func, const vector<NamedParamSpec> &params) {
	DeclareNamedParameters(func, params);
}

//! The value bound to `name`, or nullptr when the call left it out.
inline const Value *FindNamedParameter(const named_parameter_map_t &params, const char *name) {
	auto entry = params.find(name);
	if (entry == params.end()) {
		return nullptr;
	}
	return &entry->second;
}

#endif

//! The value bound to `name`. Throws like map::at did, so the call sites that
//! already knew the parameter was present keep their shape.
template <class MAP>
inline const Value &NamedParamAt(const MAP &params, const char *name) {
	auto value = FindNamedParameter(params, name);
	if (!value) {
		throw InvalidInputException("named parameter '%s' was not provided", name);
	}
	return *value;
}

//! Whether the call provided `name` at all.
template <class MAP>
inline bool HasNamedParam(const MAP &params, const char *name) {
	return FindNamedParameter(params, name) != nullptr;
}

//! Whether a REGISTERED table function declares a named parameter `name`.
//! v2.0 holds declared params as typed kwargs on the signature; v1.5 in the flat
//! named_parameters map. (plinking addition beyond the fleet copy — used by
//! plinking_duck_capabilities to introspect the live catalog. Kept here, next to
//! the registration shim, so the two stay in sync.)
#ifdef SITTING_DUCK_HAS_TYPED_KWARGS
inline bool TableFunctionDeclaresParam(const TableFunction &func, const string &name) {
	auto kwargs = func.GetSignature().GetTypedKwargs();
	return kwargs && kwargs->Find(Identifier(name)) != nullptr;
}
#else
inline bool TableFunctionDeclaresParam(const TableFunction &func, const string &name) {
	return func.named_parameters.find(name) != func.named_parameters.end();
}
#endif

//! A CatalogEntry's overload list changed element type: vector<TableFunction> on
//! v1.5, vector<shared_ptr<const TableFunction>> on v2.0. Overloads (not a
//! template) so each line picks the one that matches and the other is simply
//! unused; both named types are valid on both lines.
inline const TableFunction &CompatDerefTableFunction(const TableFunction &func) {
	return func;
}
inline const TableFunction &CompatDerefTableFunction(const shared_ptr<const TableFunction> &func) {
	return *func;
}

} // namespace duckdb
