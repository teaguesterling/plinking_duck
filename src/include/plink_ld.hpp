#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_ld table function with DuckDB.
void RegisterPlinkLd(ExtensionLoader &loader);

} // namespace duckdb
