#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! Register read_pfile_catalog: parse + verify + prune a .pfile.json catalog,
//! returning one row per surviving shard. plinking_duck owns pfile_catalog READING
//! (it already owns writing via the sharded pfile writer).
void RegisterPfileCatalogReader(ExtensionLoader &loader);

} // namespace duckdb
