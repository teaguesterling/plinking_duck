#pragma once

#include "duckdb.hpp"

namespace duckdb {

//! Register the plink_glm table function with DuckDB.
void RegisterPlinkGlm(ExtensionLoader &loader);

} // namespace duckdb
