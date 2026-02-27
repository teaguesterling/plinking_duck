#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_hardy table function with DuckDB.
void RegisterPlinkHardy(ExtensionLoader &loader);

} // namespace duckdb
