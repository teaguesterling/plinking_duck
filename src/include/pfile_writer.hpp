#pragma once

#include "duckdb.hpp"

namespace duckdb {

class ExtensionLoader;

//! Register the pfile family of COPY writers.
//!
//! P9 Phase A, first increment: `COPY (query) TO 'path' (FORMAT pvar)` — a
//! text `.pvar` writer. `.psam`, `.pgen`, and the bundled `pfile` fan-out are
//! later increments and are NOT registered here yet.
void RegisterPfileWriter(ExtensionLoader &loader);

} // namespace duckdb
