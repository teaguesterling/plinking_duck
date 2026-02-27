#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_score table function with DuckDB.
void RegisterPlinkScore(ExtensionLoader &loader);

} // namespace duckdb
