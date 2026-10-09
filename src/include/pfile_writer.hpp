#pragma once

#include "duckdb.hpp"

namespace duckdb {

class ExtensionLoader;

//! Register the pfile family of COPY writers.
//!
//! Registers:
//!   * `COPY (query) TO 'path' (FORMAT pvar)` — text `.pvar` writer (P9 Phase A).
//!   * `COPY (query) TO 'path' (FORMAT pgen)` — binary `.pgen` writer, biallelic
//!     HARDCALLS only (P9 Phase B). Phased/dosage/multiallelic are rejected.
//!     Accepts LIST/ARRAY or STRUCT genotypes (STRUCT field order = sample order).
//!   * `COPY (query) TO 'path' (FORMAT psam)` — text `.psam` writer, IID-only
//!     (P9 Phase C). Identity from STRUCT genotypes field names, a `samples`
//!     list, or a per-row IID column.
//!   * `COPY (query) TO 'prefix' (FORMAT pfile)` — the bundle (P9 Phase C): fans
//!     out `prefix.pvar` + `prefix.pgen` + `prefix.psam` from one query. Local-only
//!     (`.pgen` involved); sample identity via STRUCT field names or `samples`.
void RegisterPfileWriter(ExtensionLoader &loader);

} // namespace duckdb
