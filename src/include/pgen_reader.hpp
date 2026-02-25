#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the read_pgen table function with DuckDB.
void RegisterPgenReader(ExtensionLoader &loader);

} // namespace duckdb
