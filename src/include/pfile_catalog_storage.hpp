#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! Register the `pfile_catalog` storage extension, so a matrix catalog attaches as a
//! read-only database: `ATTACH 'cohort.pfile.json' AS c (TYPE pfile_catalog)`. The
//! attached catalog is a native DuckCatalog seeded with lazy views over the shipped
//! read_pfile_catalog / read_pfile readers (P13-003). plinking_duck already owns the
//! catalog format (writing + reading); this adds only the attach ergonomics.
void RegisterPfileCatalogStorage(ExtensionLoader &loader);

} // namespace duckdb
