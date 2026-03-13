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
#include "plink_glm.hpp"
#include "duckdb.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/config.hpp"

namespace duckdb {

// ---------------------------------------------------------------------------
// Extension config option callbacks
// ---------------------------------------------------------------------------

static void SetPlinkingMaxMatrixElements(ClientContext &, SetScope, Value &parameter) {
	// Validation only — the value is stored by DuckDB's config system
	// and retrieved via TryGetCurrentSetting at bind time.
	auto val = parameter.GetValue<int64_t>();
	if (val <= 0) {
		throw InvalidInputException("plinking_max_matrix_elements must be positive");
	}
}

void PlinkingDuckExtension::Load(ExtensionLoader &loader) {
	// Register config options
	auto &db = loader.GetDatabaseInstance();
	auto &config = DBConfig::GetConfig(db);

	config.AddExtensionOption("plinking_max_matrix_elements",
	                          "Maximum genotype matrix elements for orient := 'sample' pre-read "
	                          "(variants x samples). Default 16 billion (~16 GB of int8).",
	                          LogicalType::BIGINT, Value::BIGINT(16LL * 1024 * 1024 * 1024),
	                          SetPlinkingMaxMatrixElements);

	// Register table functions
	RegisterPvarReader(loader);
	RegisterPsamReader(loader);
	RegisterPgenReader(loader);
	RegisterPfileReader(loader);
	RegisterPlinkFreq(loader);
	RegisterPlinkHardy(loader);
	RegisterPlinkMissing(loader);
	RegisterPlinkLd(loader);
	RegisterPlinkScore(loader);
	RegisterPlinkGlm(loader);
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
