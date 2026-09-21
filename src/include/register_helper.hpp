#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

namespace duckdb {

inline void RegisterScalarWithDesc(ExtensionLoader &loader, ScalarFunction fn, const vector<string> &params = {},
                                   const string &desc_str = "", const vector<string> &examples = {}) {
	CreateScalarFunctionInfo info(std::move(fn));
	info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	FunctionDescription desc;
	desc.parameter_names = params;
	desc.description = desc_str;
	desc.examples = examples;
	desc.categories = {"plinking_duck"};
	info.descriptions.push_back(desc);
	loader.RegisterFunction(std::move(info));
}

inline void RegisterScalarSetWithDesc(ExtensionLoader &loader, ScalarFunctionSet fn_set,
                                      const vector<string> &params = {}, const string &desc_str = "",
                                      const vector<string> &examples = {}) {
	CreateScalarFunctionInfo info(std::move(fn_set));
	info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	FunctionDescription desc;
	desc.parameter_names = params;
	desc.description = desc_str;
	desc.examples = examples;
	desc.categories = {"plinking_duck"};
	info.descriptions.push_back(desc);
	loader.RegisterFunction(std::move(info));
}

inline void RegisterTableWithDesc(ExtensionLoader &loader, TableFunction fn, const vector<string> &params = {},
                                  const string &desc_str = "", const vector<string> &examples = {}) {
	CreateTableFunctionInfo info(std::move(fn));
	info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	FunctionDescription desc;
	desc.parameter_names = params;
	desc.description = desc_str;
	desc.examples = examples;
	desc.categories = {"plinking_duck"};
	info.descriptions.push_back(desc);
	loader.RegisterFunction(std::move(info));
}

inline void RegisterTableSetWithDesc(ExtensionLoader &loader, TableFunctionSet fn_set,
                                     const vector<string> &params = {}, const string &desc_str = "",
                                     const vector<string> &examples = {}) {
	CreateTableFunctionInfo info(std::move(fn_set));
	info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	FunctionDescription desc;
	desc.parameter_names = params;
	desc.description = desc_str;
	desc.examples = examples;
	desc.categories = {"plinking_duck"};
	info.descriptions.push_back(desc);
	loader.RegisterFunction(std::move(info));
}

} // namespace duckdb
