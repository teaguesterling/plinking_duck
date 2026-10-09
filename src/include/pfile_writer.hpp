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
//!
//! `.psam` and the bundled `pfile` fan-out are later increments and are NOT
//! registered here yet.
void RegisterPfileWriter(ExtensionLoader &loader);

} // namespace duckdb
