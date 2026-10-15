#pragma once

#include "duckdb.hpp"
#include "duckdb_compat.hpp"

#include "duckdb/planner/table_filter.hpp"
#include "duckdb/planner/filter/constant_filter.hpp"
#include "duckdb/planner/filter/in_filter.hpp"
#include "duckdb/planner/filter/conjunction_filter.hpp"
#include "duckdb/planner/filter/optional_filter.hpp"
#include "duckdb/planner/expression/bound_conjunction_expression.hpp"

#include <type_traits>

// Planner-side compatibility shims, kept OUT of duckdb_compat.hpp on purpose.
//
// duckdb_compat.hpp is the header the whole duck_block fleet shares and merges by
// name; it must stay cheap enough to include everywhere. These shims need the
// planner/filter-pushdown headers, which only the filter-pushdown reader wants, so
// they live here instead. If another extension in the fleet grows filter pushdown,
// this is the piece to lift.
//
// LINEAGE: lifted verbatim from teaguesterling/duckdb_zim's
// src/include/duckdb_compat_planner.hpp (the author pre-designated it as the fleet
// donor). plinking_duck ADDED, beyond zim's version: CompatConjunctionAndFilter,
// CompatConjunctionOrFilter, CompatOptionalFilter and their ::TYPE-derived
// enumerators (plinking's read_pfile pushdown recognizes AND/OR/OPTIONAL shapes;
// zim only needed CONSTANT/IN). If zim's copy grows these, de-duplicate toward it.
//
// Same rule as the shared header: FEATURE DETECTION, and probe for the API that
// exists only on v2.0 -- several of the v1.5 spellings were deprecated rather
// than deleted, so probing for the OLD one is true on both lines and silently
// never takes the new branch.

// v2.0 moved TableFilterSet out of table_filter.hpp into its own header AND
// renamed the concrete legacy filter classes in the same change. This one probe
// therefore gates both, which is normally the thing to avoid -- but there is no
// SFINAE probe for the mere existence of a type NAME, and if the two ever did
// come apart the result here is a compile error naming the missing symbol, not a
// silently-wrong branch. A loud failure is an acceptable price; a quiet one is
// not.
#if __has_include("duckdb/planner/table_filter_set.hpp")
#define DUCKDB_HAS_LEGACY_TABLE_FILTERS 1
#include "duckdb/planner/table_filter_set.hpp"
#endif

namespace duckdb {

// --- concrete pushed-down filter kinds -----------------------------------------
// v1.5: ConstantFilter / InFilter, TableFilterType::CONSTANT_COMPARISON / IN_FILTER
// v2.0: the same classes renamed LegacyConstantFilter / LegacyInFilter and the
//       enumerators prefixed LEGACY_, because v2.0's planner represents filters as
//       a general ExpressionFilter and keeps these only for deserialization and
//       expression conversion.
//
// IMPORTANT, and the reason this is a rename and not a port: on v2.0 a running
// scan is handed an ExpressionFilter, so code that matches on these legacy kinds
// simply never matches. In a best-effort pruner that costs the pushdown
// OPTIMIZATION and nothing else -- correctness must live in a residual applicator
// that re-applies every filter through the virtual TableFilter::ToExpression,
// which ExpressionFilter implements. So the v2.0 build returns the same rows,
// reading more than it strictly must. Teaching the pushdown to unwrap an
// ExpressionFilter is a real feature, not a compatibility shim, and is left out.
#ifdef DUCKDB_HAS_LEGACY_TABLE_FILTERS
using CompatConstantFilter = LegacyConstantFilter;
using CompatInFilter = LegacyInFilter;
using CompatConjunctionAndFilter = LegacyConjunctionAndFilter; // plinking addition
using CompatConjunctionOrFilter = LegacyConjunctionOrFilter;   // plinking addition
using CompatOptionalFilter = LegacyOptionalFilter;             // plinking addition
#else
using CompatConstantFilter = ConstantFilter;
using CompatInFilter = InFilter;
using CompatConjunctionAndFilter = ConjunctionAndFilter; // plinking addition
using CompatConjunctionOrFilter = ConjunctionOrFilter;   // plinking addition
using CompatOptionalFilter = OptionalFilter;             // plinking addition
#endif

// Taken from the class itself rather than spelled out, so the enumerator name is
// named in exactly one place (the alias above) and cannot drift from it.
static constexpr const TableFilterType CompatConstantFilterType = CompatConstantFilter::TYPE;
static constexpr const TableFilterType CompatInFilterType = CompatInFilter::TYPE;
static constexpr const TableFilterType CompatConjunctionAndFilterType = CompatConjunctionAndFilter::TYPE; // plinking
static constexpr const TableFilterType CompatConjunctionOrFilterType = CompatConjunctionOrFilter::TYPE;   // plinking
static constexpr const TableFilterType CompatOptionalFilterType = CompatOptionalFilter::TYPE;             // plinking

// --- iterating a TableFilterSet -------------------------------------------------
// v1.5: `filters` is a public map<idx_t, unique_ptr<TableFilter>>; the set has no
//       begin()/end().
// v2.0: `filters` is private and keyed by ProjectionIndex; the set exposes
//       begin()/end() yielding entries with GetIndex() / Filter().
//
// Flattening to a vector of (projection index, filter) is the one shape both can
// produce. The sets are tiny -- one entry per filtered column -- so materialising
// them costs nothing, and it keeps the call site a plain range-for either way.
struct CompatTableFilterEntry {
	idx_t index;
	reference<TableFilter> filter;
};

template <class T, class = void>
struct CompatFilterSetHasBegin : std::false_type {};
template <class T>
struct CompatFilterSetHasBegin<T, decltype(void(std::declval<T &>().begin()))> : std::true_type {};

template <class SET>
inline vector<CompatTableFilterEntry> CompatTableFiltersImpl(SET &set, std::true_type) {
	vector<CompatTableFilterEntry> result;
	for (auto &entry : set) {
		result.push_back(CompatTableFilterEntry {entry.GetIndex().GetIndex(), entry.Filter()});
	}
	return result;
}
template <class SET>
inline vector<CompatTableFilterEntry> CompatTableFiltersImpl(SET &set, std::false_type) {
	vector<CompatTableFilterEntry> result;
	for (auto &entry : set.filters) {
		result.push_back(CompatTableFilterEntry {entry.first, *entry.second});
	}
	return result;
}
inline vector<CompatTableFilterEntry> CompatTableFilters(TableFilterSet &set) {
	return CompatTableFiltersImpl(set, CompatFilterSetHasBegin<TableFilterSet>());
}

// v2.0 added a SECOND collection to the set: filters that span several columns and
// therefore have no single column index. They are not reachable through begin()/end()
// above, so a caller that walks only the per-column filters would not see them --
// and DuckDB deletes a fully-pushed filter from the plan, so not seeing one means
// not applying it. v1.5 has no such collection, hence the false branch.
//
// This exists so the caller can REFUSE rather than silently drop; there is no
// version of "ignore it" that is correct.
template <class T, class = void>
struct CompatFilterSetHasMultiColumn : std::false_type {};
template <class T>
struct CompatFilterSetHasMultiColumn<T, decltype(void(std::declval<const T &>().HasMultiColumnFilters()))>
    : std::true_type {};

template <class SET>
inline bool CompatHasMultiColumnFiltersImpl(const SET &set, std::true_type) {
	return set.HasMultiColumnFilters();
}
template <class SET>
inline bool CompatHasMultiColumnFiltersImpl(const SET &, std::false_type) {
	return false;
}
inline bool CompatHasMultiColumnFilters(const TableFilterSet &set) {
	return CompatHasMultiColumnFiltersImpl(set, CompatFilterSetHasMultiColumn<TableFilterSet>());
}

// --- per-column emptiness (const, no materialization) ---------------------------
// plinking addition (not in zim's donor): a const-correct "are there per-column
// filters?" that gates the prune/residual work without materializing the vector.
// v2.0 exposes HasFilters(); v1.5 has the public `filters` map.
template <class T, class = void>
struct CompatFilterSetHasHasFilters : std::false_type {};
template <class T>
struct CompatFilterSetHasHasFilters<T, decltype(void(std::declval<const T &>().HasFilters()))> : std::true_type {};

template <class SET>
inline bool CompatTableFilterSetHasFiltersImpl(const SET &set, std::true_type) {
	return set.HasFilters();
}
template <class SET>
inline bool CompatTableFilterSetHasFiltersImpl(const SET &set, std::false_type) {
	return !set.filters.empty();
}
inline bool CompatTableFilterSetHasFilters(const TableFilterSet &set) {
	return CompatTableFilterSetHasFiltersImpl(set, CompatFilterSetHasHasFilters<TableFilterSet>());
}

// --- BoundConjunctionExpression children ----------------------------------------
// v2.0 made `children` private and added GetChildren() / GetChildrenMutable().
// Probed on the MUTABLE accessor: that is the one this code needs, and it is the
// one that exists only on v2.0.
template <class T, class = void>
struct CompatHasGetChildrenMutable : std::false_type {};
template <class T>
struct CompatHasGetChildrenMutable<T, decltype(void(std::declval<T &>().GetChildrenMutable()))> : std::true_type {};

template <class EXPR>
inline vector<unique_ptr<Expression>> &CompatConjunctionChildrenImpl(EXPR &expr, std::true_type) {
	return expr.GetChildrenMutable();
}
template <class EXPR>
inline vector<unique_ptr<Expression>> &CompatConjunctionChildrenImpl(EXPR &expr, std::false_type) {
	return expr.children;
}
inline vector<unique_ptr<Expression>> &CompatConjunctionChildren(BoundConjunctionExpression &expr) {
	return CompatConjunctionChildrenImpl(expr, CompatHasGetChildrenMutable<BoundConjunctionExpression>());
}

} // namespace duckdb
