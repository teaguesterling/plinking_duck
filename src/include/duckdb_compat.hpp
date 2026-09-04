#pragma once

#include "duckdb.hpp"
#include <type_traits>
#include <utility>

// duckdb_compat.hpp — fleet-standard cross-version shim for DuckDB extensions.
//
// Pattern established by @bendrucker in teaguesterling/duckdb_webbed#76 (May 2026):
// detect the new API via __has_include of headers that moved in the same DuckDB
// refactor ([duckdb/duckdb#22377](https://github.com/duckdb/duckdb/pull/22377) —
// "mandatory per-vector size tracking" landed alongside the vector-buffer header
// reshuffle), then dispatch via a single #ifdef block.
//
// Cross-version coverage:
//   - duckdb v1.4.x / v1.5.x: old API everywhere
//   - duckdb main / v1.6.x:   new API everywhere
//
// See teaguesterling/duckdb_markdown's docs/DUCKDB_API_MIGRATION.md for the
// long-form rationale + upgrade checklist for other extensions.

#if __has_include("duckdb/common/vector/list_vector.hpp")
#define DUCKDB_HAS_NEW_VECTOR_HEADERS 1
#include "duckdb/common/vector/list_vector.hpp"
#include "duckdb/common/vector/struct_vector.hpp"
#endif

// duckdb::Identifier replaced std::string as the name type in table-function and
// COPY bind signatures (the DuckDB v2.0 / v2.0-cyanoptera line). Identifier
// compares case-insensitively, and construction from a RUNTIME string is
// explicit by design -- promoting a string to an identifier is meant to be a
// deliberate act at the call site -- so a boundary helper is needed rather than
// an implicit conversion.
//
// PROBED SEPARATELY from the vector-header change above. Each API change gets
// its own probe: tying several to one macro silently picks the wrong branch if
// they ever land in different releases (backport, revert, unexpected branch).
#if __has_include("duckdb/common/identifier.hpp")
#define DUCKDB_HAS_IDENTIFIER 1
#include "duckdb/common/identifier.hpp"
#endif

namespace duckdb {

#ifdef DUCKDB_HAS_NEW_VECTOR_HEADERS

// --- Output chunk finalization ---
// DuckDB main mandates per-vector Size() tracking; DataChunk::SetCardinality only
// updates chunk.count. SetChildCardinality additionally calls FlatVector::SetSize
// on every column so query operators reading vec.Size() see the right value.
// Without this, VariadicExecutor (and similar) reports:
//   "Mismatch in input vector sizes ... expected 0 rows but got N"
inline void CompatSetOutputCardinality(DataChunk &chunk, idx_t count) {
	chunk.SetChildCardinality(count);
}

#else // Old API (v1.4.x / v1.5.x)

inline void CompatSetOutputCardinality(DataChunk &chunk, idx_t count) {
	chunk.SetCardinality(count);
}

#endif

// --- bind-signature name type -------------------------------------------------
// Used wherever a bind callback receives or fills a vector of column names
// (table_function_bind_t, copy_to_bind_t).
//
// String literals need no help: Identifier's const char* constructor is
// IMPLICIT (a literal is an identifier by intent) while its std::string
// constructor is EXPLICIT. So `names = {"CHROM", "POS"}` and
// `names.push_back("genotypes")` compile unchanged on both lines, and only the
// signatures plus the few RUNTIME-string boundaries need the helpers below.
//
// CompatNameStr/CompatMakeName are overloaded on both string and Identifier
// where both types exist, so a call site stays correct whichever type the
// surrounding DuckDB API hands it (e.g. QueryResult::names).
#ifdef DUCKDB_HAS_IDENTIFIER
using CompatName = Identifier;
inline string CompatNameStr(const Identifier &id) {
	return id.GetIdentifierName();
}
inline string CompatNameStr(const string &name) {
	return name;
}
inline Identifier CompatMakeName(string name) {
	return Identifier(std::move(name));
}
inline Identifier CompatMakeName(const Identifier &name) {
	return name;
}
#else
using CompatName = string;
inline string CompatNameStr(const string &name) {
	return name;
}
inline string CompatMakeName(string name) {
	return name;
}
#endif

// Bulk conversion of a name vector to plain strings.
//
// Not only bind parameters moved to Identifier: `QueryResult::names` is
// `vector<Identifier>` on v2.0 too (this extension reads schemas back out of
// conn.Query() results for the non-native pvar/psam/pfile source paths). This
// template accepts vector<string> or vector<Identifier>, so a call site that
// only wants text stays correct on both lines.
template <class NAMES>
inline vector<string> CompatNameStrings(const NAMES &names) {
	vector<string> result;
	result.reserve(names.size());
	for (auto &name : names) {
		result.push_back(CompatNameStr(name));
	}
	return result;
}

// --- StructVector child entries -------------------------------------------------
// v1.5: StructVector::GetEntries(vec) -> vector<unique_ptr<Vector>> &
// v2.0: StructVector::GetEntries(vec) -> vector<Vector> &
// so the idiomatic `*entries[i]` stops compiling on v2.0 (you would be
// dereferencing a Vector). Resolved by OVERLOADING on the element type rather
// than by a probe: the overload set is unambiguous on both lines, and the
// unused overload is simply never selected.
inline Vector &CompatStructEntry(Vector &entry) {
	return entry;
}
inline Vector &CompatStructEntry(unique_ptr<Vector> &entry) {
	return *entry;
}

// --- LogicalType alias ---------------------------------------------------------
// v1.5: void SetAlias(string)                -- mutates in place
// v2.0: LogicalType WithAlias(string) const  -- returns a copy, never mutating a
//       type whose type-info may be shared. SetAlias is REMOVED, not deprecated.
//
// Detected by PROBING for the member rather than by the Identifier macro above,
// because these are independent changes. `if constexpr` discards the untaken
// branch only inside a template, hence the template parameter.
template <class T, class = void>
struct CompatHasWithAlias : std::false_type {};
template <class T>
struct CompatHasWithAlias<T, decltype(void(std::declval<const T &>().WithAlias(string())))> : std::true_type {};

template <class TYPE = LogicalType>
inline LogicalType CompatWithAlias(TYPE type, string alias) {
	if constexpr (CompatHasWithAlias<TYPE>::value) {
		return type.WithAlias(std::move(alias));
	} else {
		type.SetAlias(std::move(alias));
		return type;
	}
}

// --- Vector::ToUnifiedFormat ---------------------------------------------------
// v2.0 dropped the count parameter (the vector now tracks its own size). Probed
// the same way.
template <class T, class = void>
struct CompatToUnifiedTakesCount : std::false_type {};
template <class T>
struct CompatToUnifiedTakesCount<T, decltype(void(std::declval<T &>().ToUnifiedFormat(
                                        idx_t(0), std::declval<UnifiedVectorFormat &>())))> : std::true_type {};

template <class VEC = Vector>
inline void CompatToUnifiedFormat(VEC &vec, idx_t count, UnifiedVectorFormat &data) {
	if constexpr (CompatToUnifiedTakesCount<VEC>::value) {
		vec.ToUnifiedFormat(count, data);
	} else {
		vec.ToUnifiedFormat(data);
	}
}

// --- FlatVector mutable data ---------------------------------------------------
// v1.5: FlatVector::GetData<T>(vec)         returns T*
// v2.0: FlatVector::GetData<T>(vec)         returns const T*
//       FlatVector::GetDataMutable<T>(vec)  returns T*
// Writing through the v2.0 read accessor is a compile error, which is the point
// of the split -- so every WRITE path must ask for mutability explicitly. Read
// paths keep using FlatVector::GetData directly.
template <class T, class = void>
struct CompatHasFlatGetDataMutable : std::false_type {};
template <class T>
struct CompatHasFlatGetDataMutable<T, decltype(void(T::template GetDataMutable<bool>(std::declval<Vector &>())))>
    : std::true_type {};

template <class VALUE, class FV = FlatVector>
inline VALUE *CompatFlatDataMutable(Vector &vec) {
	if constexpr (CompatHasFlatGetDataMutable<FV>::value) {
		return FV::template GetDataMutable<VALUE>(vec);
	} else {
		return FV::template GetData<VALUE>(vec);
	}
}

// --- FlatVector mutable validity mask -------------------------------------------
// The same const split applies to the validity mask:
// v1.5: FlatVector::Validity(vec)          returns ValidityMask&
// v2.0: FlatVector::Validity(vec)          returns const ValidityMask&
//       FlatVector::ValidityMutable(vec)   returns ValidityMask&
// so SetInvalid()/SetValid() through the read accessor stops compiling.
// Probed independently of GetDataMutable above.
template <class T, class = void>
struct CompatHasFlatValidityMutable : std::false_type {};
template <class T>
struct CompatHasFlatValidityMutable<T, decltype(void(T::ValidityMutable(std::declval<Vector &>())))> : std::true_type {
};

template <class FV = FlatVector>
inline ValidityMask &CompatFlatValidityMutable(Vector &vec) {
	if constexpr (CompatHasFlatValidityMutable<FV>::value) {
		return FV::ValidityMutable(vec);
	} else {
		return FV::Validity(vec);
	}
}

} // namespace duckdb
