#include "pvar_reader.hpp"

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <limits>

namespace duckdb {

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

//! Map known .pvar column names to their DuckDB types.
//! POS is always INTEGER, QUAL is FLOAT, CM is DOUBLE.
//! All other columns (CHROM, ID, REF, ALT, FILTER, INFO, and any
//! unrecognized names) are VARCHAR.
static LogicalType PvarColumnType(const string &name) {
	if (name == "POS") {
		return LogicalType::INTEGER;
	}
	if (name == "QUAL") {
		return LogicalType::FLOAT;
	}
	if (name == "CM") {
		return LogicalType::DOUBLE;
	}
	return LogicalType::VARCHAR;
}

PvarHeaderInfo ParsePvarHeader(const string &file_path) {
	std::ifstream file(file_path);
	if (!file.is_open()) {
		throw IOException("read_pvar: could not open file '%s'", file_path);
	}

	PvarHeaderInfo info;
	info.skip_lines = 0;
	info.is_bim = false;
	string line;
	bool found_header_or_data = false;

	while (std::getline(file, line)) {
		// Handle Windows line endings
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		// Skip empty lines
		if (line.empty()) {
			info.skip_lines++;
			continue;
		}
		// Skip ## comment/meta lines (e.g. ##fileformat=PVARv1.0)
		if (line.size() >= 2 && line[0] == '#' && line[1] == '#') {
			info.skip_lines++;
			continue;
		}
		found_header_or_data = true;
		break;
	}

	if (!found_header_or_data) {
		throw InvalidInputException("read_pvar: file '%s' is empty or contains no header/data", file_path);
	}

	// Format detection: #CHROM header indicates .pvar, otherwise .bim
	if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
		// .pvar format: parse column names from the header line
		info.is_bim = false;
		info.skip_lines++; // count the header line itself

		// Strip the leading '#' so "#CHROM" becomes "CHROM"
		string header = line.substr(1);
		size_t start = 0;
		size_t pos;
		while ((pos = header.find('\t', start)) != string::npos) {
			string col_name = header.substr(start, pos - start);
			info.column_names.push_back(col_name);
			info.column_types.push_back(PvarColumnType(col_name));
			start = pos + 1;
		}
		// Last column (no trailing tab)
		string col_name = header.substr(start);
		info.column_names.push_back(col_name);
		info.column_types.push_back(PvarColumnType(col_name));
	} else {
		// No #CHROM header: legacy .bim format with 6 fixed columns.
		//
		// .bim file column order: CHROM(0), ID(1), CM(2), POS(3), ALT(4), REF(5)
		// Normalized output order: CHROM(0), POS(1), ID(2), REF(3), ALT(4), CM(5)
		//
		// The reader normalizes .bim output to match .pvar column ordering
		// so that downstream queries work identically on both formats.
		info.is_bim = true;
		info.skip_lines = 0; // no header to skip; data starts at line 1
		info.column_names = {"CHROM", "POS", "ID", "REF", "ALT", "CM"};
		info.column_types = {
		    LogicalType::VARCHAR,  // CHROM
		    LogicalType::INTEGER,  // POS
		    LogicalType::VARCHAR,  // ID
		    LogicalType::VARCHAR,  // REF
		    LogicalType::VARCHAR,  // ALT
		    LogicalType::DOUBLE    // CM
		};
	}

	return info;
}

// ---------------------------------------------------------------------------
// Table function data structures
// ---------------------------------------------------------------------------

struct PvarBindData : public TableFunctionData {
	string file_path;
	PvarHeaderInfo header_info;
};

struct PvarGlobalState : public GlobalTableFunctionState {
	std::ifstream file;
	bool finished = false;
	vector<column_t> column_ids;

	idx_t MaxThreads() const override {
		return 1; // Sequential text file reading
	}
};

struct PvarLocalState : public LocalTableFunctionState {};

// ---------------------------------------------------------------------------
// .bim column order mapping
// ---------------------------------------------------------------------------

// Maps normalized output column index to .bim file field index.
// Output: CHROM(0) POS(1) ID(2) REF(3) ALT(4) CM(5)
// File:   CHROM(0) ID(1)  CM(2) POS(3) ALT(4) REF(5)
static constexpr idx_t BIM_OUTPUT_TO_FILE[] = {0, 3, 1, 5, 4, 2};

// ---------------------------------------------------------------------------
// Table function callbacks
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PvarBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PvarBindData>();
	bind_data->file_path = input.inputs[0].GetValue<string>();
	bind_data->header_info = ParsePvarHeader(bind_data->file_path);

	names = bind_data->header_info.column_names;
	return_types = bind_data->header_info.column_types;

	return std::move(bind_data);
}

static unique_ptr<GlobalTableFunctionState> PvarInitGlobal(ClientContext &context,
                                                           TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PvarBindData>();
	auto state = make_uniq<PvarGlobalState>();

	state->file.open(bind_data.file_path);
	if (!state->file.is_open()) {
		throw IOException("read_pvar: could not open file '%s'", bind_data.file_path);
	}

	// Skip comment and header lines to reach the first data line
	string skip;
	for (idx_t i = 0; i < bind_data.header_info.skip_lines; i++) {
		std::getline(state->file, skip);
	}

	// Store projected column indices for the scan function
	state->column_ids = input.column_ids;

	return std::move(state);
}

static unique_ptr<LocalTableFunctionState> PvarInitLocal(ExecutionContext &context,
                                                         TableFunctionInitInput &input,
                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<PvarLocalState>();
}

// ---------------------------------------------------------------------------
// Scan helpers
// ---------------------------------------------------------------------------

//! Split a line on tab characters into fields.
//! Used for .pvar format where tabs are the only valid delimiter
//! (fields like INFO can contain spaces).
static vector<string> SplitTabs(const string &line) {
	vector<string> fields;
	size_t start = 0;
	size_t pos;
	while ((pos = line.find('\t', start)) != string::npos) {
		fields.push_back(line.substr(start, pos - start));
		start = pos + 1;
	}
	fields.push_back(line.substr(start));
	return fields;
}

//! Split a line on whitespace (spaces and/or tabs) into fields.
//! Used for .bim format where PLINK 1 allows any whitespace delimiter.
//! Consecutive whitespace is treated as a single delimiter.
static vector<string> SplitWhitespace(const string &line) {
	vector<string> fields;
	size_t i = 0;
	while (i < line.size()) {
		// Skip whitespace between fields
		while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) {
			i++;
		}
		if (i >= line.size()) {
			break;
		}
		// Read field until next whitespace
		size_t start = i;
		while (i < line.size() && line[i] != ' ' && line[i] != '\t') {
			i++;
		}
		fields.push_back(line.substr(start, i - start));
	}
	return fields;
}

//! Parse a field value and write it to an output vector.
//! A single dot (".") is treated as NULL for any column type.
static void SetPvarValue(Vector &vec, idx_t row_idx, const string &field, const LogicalType &type) {
	if (field == ".") {
		FlatVector::Validity(vec).SetInvalid(row_idx);
		return;
	}

	switch (type.id()) {
	case LogicalTypeId::VARCHAR: {
		FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, field);
		break;
	}
	case LogicalTypeId::INTEGER: {
		char *end;
		errno = 0;
		long val = std::strtol(field.c_str(), &end, 10);
		if (end == field.c_str() || *end != '\0' || errno != 0) {
			throw InvalidInputException("read_pvar: invalid integer value '%s'", field);
		}
		if (val > static_cast<long>(std::numeric_limits<int32_t>::max()) ||
		    val < static_cast<long>(std::numeric_limits<int32_t>::min())) {
			throw InvalidInputException("read_pvar: integer value '%s' out of range", field);
		}
		FlatVector::GetData<int32_t>(vec)[row_idx] = static_cast<int32_t>(val);
		break;
	}
	case LogicalTypeId::FLOAT: {
		char *end;
		errno = 0;
		float val = std::strtof(field.c_str(), &end);
		if (end == field.c_str() || *end != '\0' || errno != 0) {
			throw InvalidInputException("read_pvar: invalid float value '%s'", field);
		}
		FlatVector::GetData<float>(vec)[row_idx] = val;
		break;
	}
	case LogicalTypeId::DOUBLE: {
		char *end;
		errno = 0;
		double val = std::strtod(field.c_str(), &end);
		if (end == field.c_str() || *end != '\0' || errno != 0) {
			throw InvalidInputException("read_pvar: invalid double value '%s'", field);
		}
		FlatVector::GetData<double>(vec)[row_idx] = val;
		break;
	}
	default:
		throw InternalException("read_pvar: unsupported column type");
	}
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void PvarScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PvarBindData>();
	auto &state = data_p.global_state->Cast<PvarGlobalState>();

	if (state.finished) {
		output.SetCardinality(0);
		return;
	}

	auto &header = bind_data.header_info;
	idx_t row_count = 0;
	string line;

	while (row_count < STANDARD_VECTOR_SIZE && std::getline(state.file, line)) {
		// Handle Windows line endings
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		// Skip empty lines
		if (line.empty()) {
			continue;
		}

		// .bim uses whitespace-delimited fields (spaces, tabs, or mixed);
		// .pvar uses tab-only (fields like INFO can contain spaces)
		auto fields = header.is_bim ? SplitWhitespace(line) : SplitTabs(line);

		// Validate field count
		idx_t expected = header.is_bim ? 6 : header.column_names.size();
		if (fields.size() < expected) {
			throw InvalidInputException(
			    "read_pvar: line has %llu fields, expected at least %llu in '%s'",
			    static_cast<unsigned long long>(fields.size()),
			    static_cast<unsigned long long>(expected), bind_data.file_path);
		}

		// Fill projected output columns
		for (idx_t i = 0; i < state.column_ids.size(); i++) {
			auto col_idx = state.column_ids[i];
			if (col_idx == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			// Map output column index to file field index
			idx_t file_idx = header.is_bim ? BIM_OUTPUT_TO_FILE[col_idx] : col_idx;
			SetPvarValue(output.data[i], row_count, fields[file_idx], header.column_types[col_idx]);
		}

		row_count++;
	}

	if (row_count == 0) {
		state.finished = true;
	}

	output.SetCardinality(row_count);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPvarReader(ExtensionLoader &loader) {
	TableFunction read_pvar("read_pvar", {LogicalType::VARCHAR}, PvarScan, PvarBind, PvarInitGlobal,
	                        PvarInitLocal);
	read_pvar.projection_pushdown = true;
	loader.RegisterFunction(read_pvar);
}

} // namespace duckdb
