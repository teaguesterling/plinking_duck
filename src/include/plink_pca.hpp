#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_pca table function with DuckDB.
void RegisterPlinkPca(ExtensionLoader &loader);

} // namespace duckdb
