#define DUCKDB_EXTENSION_MAIN

#include "plinking_duck_extension.hpp"
#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

// OpenSSL linked through vcpkg
#include <openssl/opensslv.h>

namespace duckdb {

inline void PlinkingDuckScalarFun(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &name_vector = args.data[0];
	UnaryExecutor::Execute<string_t, string_t>(name_vector, result, args.size(), [&](string_t name) {
		return StringVector::AddString(result, "PlinkingDuck " + name.GetString() + " 🐥");
	});
}

inline void PlinkingDuckOpenSSLVersionScalarFun(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &name_vector = args.data[0];
	UnaryExecutor::Execute<string_t, string_t>(name_vector, result, args.size(), [&](string_t name) {
		return StringVector::AddString(result, "PlinkingDuck " + name.GetString() + ", my linked OpenSSL version is " +
		                                           OPENSSL_VERSION_TEXT);
	});
}

static void LoadInternal(ExtensionLoader &loader) {
	// Register a scalar function
	auto plinking_duck_scalar_function = ScalarFunction("plinking_duck", {LogicalType::VARCHAR}, LogicalType::VARCHAR, PlinkingDuckScalarFun);
	loader.RegisterFunction(plinking_duck_scalar_function);

	// Register another scalar function
	auto plinking_duck_openssl_version_scalar_function = ScalarFunction("plinking_duck_openssl_version", {LogicalType::VARCHAR},
	                                                            LogicalType::VARCHAR, PlinkingDuckOpenSSLVersionScalarFun);
	loader.RegisterFunction(plinking_duck_openssl_version_scalar_function);
}

void PlinkingDuckExtension::Load(ExtensionLoader &loader) {
	LoadInternal(loader);
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

DUCKDB_CPP_EXTENSION_ENTRY(plinking_duck, loader) {
	duckdb::LoadInternal(loader);
}
}
