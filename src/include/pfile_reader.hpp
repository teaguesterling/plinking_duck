#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the read_pfile table function with DuckDB.
void RegisterPfileReader(ExtensionLoader &loader);

} // namespace duckdb
