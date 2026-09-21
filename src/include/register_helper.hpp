#pragma once

#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"

namespace duckdb {

inline void RegisterScalarWithDesc(ExtensionLoader &loader, ScalarFunction func, vector<string> parameter_names,
                                   string description, vector<string> examples = {}) {
	FunctionDescription desc(parameter_names, description, examples, {"plinking_duck"});
	func.descriptions.push_back(std::move(desc));
	auto info = make_uniq<CreateScalarFunctionInfo>(std::move(func));
	info->on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	loader.RegisterFunction(std::move(info));
}

inline void RegisterScalarSetWithDesc(ExtensionLoader &loader, ScalarFunctionSet set, vector<string> parameter_names,
                                      string description, vector<string> examples = {}) {
	FunctionDescription desc(parameter_names, description, examples, {"plinking_duck"});
	set.descriptions.push_back(std::move(desc));
	auto info = make_uniq<CreateScalarFunctionInfo>(std::move(set));
	info->on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	loader.RegisterFunction(std::move(info));
}

inline void RegisterTableWithDesc(ExtensionLoader &loader, TableFunction func, vector<string> parameter_names,
                                  string description, vector<string> examples = {}) {
	FunctionDescription desc(parameter_names, description, examples, {"plinking_duck"});
	func.descriptions.push_back(std::move(desc));
	auto info = make_uniq<CreateTableFunctionInfo>(std::move(func));
	info->on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	loader.RegisterFunction(std::move(info));
}

inline void RegisterTableSetWithDesc(ExtensionLoader &loader, TableFunctionSet set, vector<string> parameter_names,
                                     string description, vector<string> examples = {}) {
	FunctionDescription desc(parameter_names, description, examples, {"plinking_duck"});
	set.descriptions.push_back(std::move(desc));
	auto info = make_uniq<CreateTableFunctionInfo>(std::move(set));
	info->on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	loader.RegisterFunction(std::move(info));
}

} // namespace duckdb
