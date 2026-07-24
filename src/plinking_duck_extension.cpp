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
#include "vcf_reader.hpp"
#ifdef PLINKING_HAVE_EIGEN3
#include "plink_pca.hpp"
#endif
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

static void SetPlinkingMaxThreads(ClientContext &, SetScope, Value &parameter) {
	auto val = parameter.GetValue<int64_t>();
	if (val < 0) {
		throw InvalidInputException("plinking_max_threads must be non-negative (0 = default, >0 = cap)");
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

	config.AddExtensionOption("plinking_max_threads",
	                          "Maximum threads for parallel scan operations. "
	                          "0 = default (hardcoded cap of 16), >0 = cap at this value.",
	                          LogicalType::BIGINT, Value::BIGINT(0), SetPlinkingMaxThreads);

	config.AddExtensionOption("plinking_use_parquet_companions",
	                          "Auto-discover .pvar.parquet and .psam.parquet companion files. "
	                          "When true, parquet companions are preferred over text formats.",
	                          LogicalType::BOOLEAN, Value::BOOLEAN(true));

	config.AddExtensionOption(
	    "plinking_pgen_io",
	    "How .pgen bytes are read: 'auto' (default — remote/VFS paths via DuckDB's VFS, local via native "
	    "fopen), 'native' (always native fopen; errors on remote), 'vfs' (always via DuckDB's VFS, even "
	    "local), 'localize' (materialize a local temp copy then read natively — best for remote full scans; "
	    "always copies, even a local source).",
	    LogicalType::VARCHAR, Value("auto"));

	config.AddExtensionOption("plinking_localize_dir",
	                          "Directory (local) for temp copies made by plinking_pgen_io := 'localize'. Empty "
	                          "(default) uses DuckDB's temporary_directory, else the current directory. Created if "
	                          "absent (one level only — parent dirs must already exist). Temps are per-query and "
	                          "removed when the query's bind data is destroyed.",
	                          LogicalType::VARCHAR, Value(""));

	config.AddExtensionOption("plinking_sample_counts_sparse",
	                          "orient := 'sample' + genotypes := 'counts'|'stats': when true, use the "
	                          "sparse (pgen difflist) accumulation path — reads only the non-hom_ref "
	                          "carriers of rare variants (auto-falls-back to the dense full-decode path "
	                          "per variant). When false, always use the dense path. Toggle to A/B time "
	                          "both; both paths produce identical counts.",
	                          LogicalType::BOOLEAN, Value::BOOLEAN(false));

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
	RegisterPlinkVcfReader(loader);
#ifdef PLINKING_HAVE_EIGEN3
	RegisterPlinkPca(loader);
#endif
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
