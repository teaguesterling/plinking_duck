#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// Registers the scalar function plinking_duck_capabilities() -> VARCHAR[].
void RegisterPlinkingCapabilities(ExtensionLoader &loader);

} // namespace duckdb
