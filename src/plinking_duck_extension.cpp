#define DUCKDB_EXTENSION_MAIN

#include "plinking_duck_extension.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"
#include "pgen_reader.hpp"
#include "pfile_reader.hpp"
#include "plink_freq.hpp"
#include "plink_hardy.hpp"
#include "plink_missing.hpp"
#include "plink_ld.hpp"
#include "plink_score.hpp"
#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

void PlinkingDuckExtension::Load(ExtensionLoader &loader) {
	RegisterPvarReader(loader);
	RegisterPsamReader(loader);
	RegisterPgenReader(loader);
	RegisterPfileReader(loader);
	RegisterPlinkFreq(loader);
	RegisterPlinkHardy(loader);
	RegisterPlinkMissing(loader);
	RegisterPlinkLd(loader);
	RegisterPlinkScore(loader);
}

std::string PlinkingDuckExtension::Name() {
	return "plinking_duck";
}

std::string PlinkingDuckExtension::Version() const {
#ifdef EXT_VERSION_PLINKING_DUCK
	return EXT_VERSION_PLINKING_DUCK;
#else
	return "";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void plinking_duck_init(duckdb::DatabaseInstance &db) {
	duckdb::DuckDB db_wrapper(db);
	db_wrapper.LoadStaticExtension<duckdb::PlinkingDuckExtension>();
}

DUCKDB_EXTENSION_API const char *plinking_duck_version() {
	return duckdb::DuckDB::LibraryVersion();
}
}

#ifdef DUCKDB_BUILD_LOADABLE_EXTENSION
extern "C" {
DUCKDB_CPP_EXTENSION_ENTRY(plinking_duck, loader) { // NOLINT
	duckdb::PlinkingDuckExtension extension;
	extension.Load(loader);
}
}
#endif
