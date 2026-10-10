#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! Register read_howard_manifest: parse + verify + prune a .howard.json catalog,
//! returning one row per surviving shard. plinking_duck owns .howard READING
//! (it already owns writing via the sharded pfile writer).
void RegisterHowardReader(ExtensionLoader &loader);

} // namespace duckdb
