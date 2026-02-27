#define DUCKDB_EXTENSION_MAIN

#include "plinking_duck_extension.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"
#include "pgen_reader.hpp"
#include "pfile_reader.hpp"
#include "plink_freq.hpp"
#include "plink_score.hpp"
#include "duckdb.hpp"

namespace duckdb {

void PlinkingDuckExtension::Load(ExtensionLoader &loader) {
	RegisterPvarReader(loader);
	RegisterPsamReader(loader);
	RegisterPgenReader(loader);
	RegisterPfileReader(loader);
	RegisterPlinkFreq(loader);
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
