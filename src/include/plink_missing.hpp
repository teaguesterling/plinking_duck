#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_missing table function with DuckDB.
void RegisterPlinkMissing(ExtensionLoader &loader);

} // namespace duckdb
