#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_freq table function with DuckDB.
void RegisterPlinkFreq(ExtensionLoader &loader);

} // namespace duckdb
