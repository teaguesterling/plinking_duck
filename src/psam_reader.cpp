#include "psam_reader.hpp"
#include "plink_common.hpp"
#include "duckdb.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <string>

namespace duckdb {

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

// .fam fixed column names (6 columns, no header line)
static const vector<string> FAM_COLUMN_NAMES = {"FID", "IID", "PAT", "MAT", "SEX", "PHENO1"};
static constexpr idx_t FAM_COLUMN_COUNT = 6;

// Known columns that have special type handling
static const string COL_SEX = "SEX";

// Missing value sentinels in .psam/.fam files
static bool IsMissingValue(const string &val) {
	return val.empty() || val == "." || val == "NA" || val == "na";
}

// PAT/MAT use "0" as missing (meaning "unknown parent")
static bool IsParentMissing(const string &val) {
	return val == "0" || IsMissingValue(val);
}

// Note on .fam PHENO1: PLINK conventionally uses -9 as a missing phenotype
// sentinel, but we intentionally leave it as the string "-9" rather than
// mapping to NULL. This avoids imposing PLINK-specific conventions on what
// is otherwise a generic VARCHAR column — users can filter with
// WHERE PHENO1 != '-9' if needed.

// ---------------------------------------------------------------------------
// Line splitting utilities (use shared implementations from plink_common)
// ---------------------------------------------------------------------------

//! Split a line using the appropriate delimiter for the given format.
static vector<string> SplitLine(const string &line, PsamFormat format) {
	if (format == PsamFormat::FAM) {
		return SplitWhitespaceLine(line);
	}
	return SplitTabLine(line);
}

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

//! Core header parsing logic, operating on already-read lines.
//! Separated from file I/O so callers that already have the lines
//! (LoadSampleInfo, PsamInitGlobal) don't re-read the file.
PsamHeaderInfo ParsePsamHeaderFromLines(const vector<string> &lines, const string &path) {
	if (lines.empty()) {
		throw IOException("read_psam: file '%s' is empty", path);
	}

	auto &first_line = lines[0];
	if (first_line.empty()) {
		throw IOException("read_psam: file '%s' has an empty first line", path);
	}

	PsamHeaderInfo info;

	if (first_line[0] == '#') {
		// .psam format: header line starting with # describes columns
		// Skip the leading '#' without mutating the input
		string header_content = first_line.substr(1);
		auto fields = SplitTabLine(header_content);

		if (fields.empty()) {
			throw IOException("read_psam: file '%s' has an empty header", path);
		}

		// Detect #FID vs #IID
		if (fields[0] == "FID") {
			info.format = PsamFormat::PSAM_FID;
		} else if (fields[0] == "IID") {
			info.format = PsamFormat::PSAM_IID;
		} else {
			throw IOException("read_psam: file '%s' header must start with #FID or #IID, got '#%s'", path, fields[0]);
		}

		for (auto &name : fields) {
			info.column_names.push_back(name);
			// SEX is the only column with a non-VARCHAR type
			if (name == COL_SEX) {
				info.column_types.push_back(LogicalType::INTEGER);
			} else {
				info.column_types.push_back(LogicalType::VARCHAR);
			}
		}
	} else {
		// .fam format: no header, fixed 6-column layout
		info.format = PsamFormat::FAM;
		for (idx_t i = 0; i < FAM_COLUMN_COUNT; i++) {
			info.column_names.push_back(FAM_COLUMN_NAMES[i]);
			if (FAM_COLUMN_NAMES[i] == COL_SEX) {
				info.column_types.push_back(LogicalType::INTEGER);
			} else {
				info.column_types.push_back(LogicalType::VARCHAR);
			}
		}
	}

	return info;
}

PsamHeaderInfo ParsePsamHeader(ClientContext &context, const string &path) {
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	auto file_size = handle->GetFileSize();

	if (file_size == 0) {
		throw IOException("read_psam: file '%s' is empty", path);
	}

	// Read the first line to determine format and column schema.
	// Use a buffered read — header lines are short (typically < 1KB).
	string first_line;
	static constexpr size_t HEADER_BUF_SIZE = 4096;
	char buf[HEADER_BUF_SIZE];
	size_t total_read = 0;
	bool found_newline = false;

	while (!found_newline && total_read < file_size) {
		size_t to_read = std::min<size_t>(HEADER_BUF_SIZE, file_size - total_read);
		auto bytes = handle->Read(buf, to_read);
		if (bytes == 0) {
			break;
		}
		for (size_t i = 0; i < bytes; i++) {
			if (buf[i] == '\n') {
				found_newline = true;
				total_read += i + 1;
				break;
			}
			if (buf[i] != '\r') {
				first_line += buf[i];
			}
		}
		if (!found_newline) {
			total_read += bytes;
		}
	}

	// Parse header from the first line
	vector<string> lines = {first_line};
	auto info = ParsePsamHeaderFromLines(lines, path);

	// For .psam: data starts after the header line
	if (info.format != PsamFormat::FAM) {
		info.data_start_offset = total_read;
	}
	// For .fam: data_start_offset stays 0 (no header)

	return info;
}

// ---------------------------------------------------------------------------
// LoadSampleInfo — reusable utility for read_pgen / read_pfile
// ---------------------------------------------------------------------------

SampleInfo LoadSampleInfo(ClientContext &context, const string &path) {
	// Single file read — ParsePsamHeaderFromLines reuses the same lines
	auto lines = ReadFileLines(context, path);
	auto header = ParsePsamHeaderFromLines(lines, path);

	// Find IID and FID column indices
	idx_t iid_idx = DConstants::INVALID_INDEX;
	idx_t fid_idx = DConstants::INVALID_INDEX;
	for (idx_t i = 0; i < header.column_names.size(); i++) {
		if (header.column_names[i] == "IID") {
			iid_idx = i;
		} else if (header.column_names[i] == "FID") {
			fid_idx = i;
		}
	}

	if (iid_idx == DConstants::INVALID_INDEX) {
		throw IOException("read_psam: file '%s' has no IID column", path);
	}

	SampleInfo info;
	bool has_fid = (fid_idx != DConstants::INVALID_INDEX);

	// Data lines start at index 1 for .psam (skip header), 0 for .fam
	idx_t data_start = (header.format != PsamFormat::FAM) ? 1 : 0;

	for (idx_t i = data_start; i < lines.size(); i++) {
		auto &line = lines[i];
		if (line.empty()) {
			continue;
		}

		auto fields = SplitLine(line, header.format);
		if (iid_idx >= fields.size()) {
			throw IOException("read_psam: file '%s' line %d has %d fields, expected at least %d", path, i + 1,
			                  fields.size(), iid_idx + 1);
		}

		const auto &iid = fields[iid_idx];

		// Duplicate IIDs are invalid — downstream code (read_pgen) relies on
		// iid_to_idx being a 1:1 mapping for sample subsetting
		if (info.iid_to_idx.count(iid)) {
			throw IOException("read_psam: file '%s' line %d has duplicate IID '%s' "
			                  "(first seen at sample %d)",
			                  path, i + 1, iid, info.iid_to_idx[iid] + 1);
		}

		info.iids.push_back(iid);
		if (has_fid) {
			info.fids.push_back(fields[fid_idx]);
		}
		info.iid_to_idx[iid] = info.iids.size() - 1;
	}

	info.sample_ct = info.iids.size();
	return info;
}

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PsamBindData : public TableFunctionData {
	string file_path;
	PsamFormat format;
	vector<string> column_names;
	vector<LogicalType> column_types;

	//! Indices of columns in the file that correspond to PAT/MAT
	//! (used for special missing-value handling where "0" means unknown)
	vector<idx_t> parent_col_indices;

	//! Index of SEX column in the file (-1 if absent)
	idx_t sex_col_idx = DConstants::INVALID_INDEX;

	//! Byte offset where data lines begin (after header)
	uint64_t data_start_offset = 0;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PsamGlobalState : public TextFileGlobalState {
	//! For each output column, whether it's a PAT/MAT column
	vector<bool> is_parent_col;
};

struct PsamLocalState : public TextFileLocalState {};

// ---------------------------------------------------------------------------
// Bind function
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PsamBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto result = make_uniq<PsamBindData>();
	result->file_path = input.inputs[0].GetValue<string>();

	auto header = ParsePsamHeader(context, result->file_path);
	result->format = header.format;
	result->column_names = header.column_names;
	result->column_types = header.column_types;
	result->data_start_offset = header.data_start_offset;

	// Identify PAT/MAT columns for special missing-value handling
	for (idx_t i = 0; i < header.column_names.size(); i++) {
		if (header.column_names[i] == "PAT" || header.column_names[i] == "MAT") {
			result->parent_col_indices.push_back(i);
		}
		if (header.column_names[i] == COL_SEX) {
			result->sex_col_idx = i;
		}
	}

	// Populate output schema
	for (idx_t i = 0; i < header.column_names.size(); i++) {
		names.push_back(header.column_names[i]);
		return_types.push_back(header.column_types[i]);
	}

	return std::move(result);
}

// ---------------------------------------------------------------------------
// Init functions
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PsamInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PsamBindData>();
	auto state = make_uniq<PsamGlobalState>();

	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(bind_data.file_path, FileFlags::FILE_FLAGS_READ);

	state->file_path = bind_data.file_path;
	state->file_size = handle->GetFileSize();
	state->data_start_offset = bind_data.data_start_offset;
	state->next_chunk_offset.store(state->data_start_offset);

	// Store projected column IDs and pre-compute parent column flags
	state->column_ids = input.column_ids;
	for (idx_t out_col = 0; out_col < state->column_ids.size(); out_col++) {
		auto file_col = state->column_ids[out_col];
		bool is_parent = false;
		for (auto &parent_idx : bind_data.parent_col_indices) {
			if (file_col == parent_idx) {
				is_parent = true;
				break;
			}
		}
		state->is_parent_col.push_back(is_parent);
	}

	return std::move(state);
}

static unique_ptr<LocalTableFunctionState> PsamInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                         GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<PsamGlobalState>();
	auto state = make_uniq<PsamLocalState>();
	state->Init(context.client, gstate.file_path);
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void PsamScan(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &bind_data = input.bind_data->Cast<PsamBindData>();
	auto &gstate = input.global_state->Cast<PsamGlobalState>();
	auto &lstate = input.local_state->Cast<PsamLocalState>();

	auto &column_ids = gstate.column_ids;
	auto &is_parent_col = gstate.is_parent_col;
	idx_t expected_cols = bind_data.column_names.size();
	idx_t row_count = 0;
	string line;

	while (row_count < STANDARD_VECTOR_SIZE) {
		// Claim a new chunk if current one is exhausted
		if (lstate.chunk_exhausted) {
			if (!lstate.ClaimChunk(gstate)) {
				break;
			}
		}

		if (!lstate.ReadChunkLine(line, gstate.file_size)) {
			continue;
		}

		if (line.empty()) {
			continue;
		}

		auto fields = SplitLine(line, bind_data.format);

		if (fields.size() != expected_cols) {
			throw IOException("read_psam: line has %llu fields, expected %llu in '%s'",
			                  static_cast<unsigned long long>(fields.size()),
			                  static_cast<unsigned long long>(expected_cols), bind_data.file_path);
		}

		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];

			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];
			const string &val = fields[file_col];

			if (file_col == bind_data.sex_col_idx) {
				if (IsMissingValue(val)) {
					FlatVector::SetNull(vec, row_count, true);
				} else {
					try {
						int32_t sex_val = std::stoi(val);
						if (sex_val == 0) {
							FlatVector::SetNull(vec, row_count, true);
						} else {
							FlatVector::GetData<int32_t>(vec)[row_count] = sex_val;
						}
					} catch (const std::invalid_argument &) {
						FlatVector::SetNull(vec, row_count, true);
					} catch (const std::out_of_range &) {
						FlatVector::SetNull(vec, row_count, true);
					}
				}
			} else if (is_parent_col[out_col]) {
				if (IsParentMissing(val)) {
					FlatVector::SetNull(vec, row_count, true);
				} else {
					FlatVector::GetData<string_t>(vec)[row_count] = StringVector::AddString(vec, val);
				}
			} else {
				if (IsMissingValue(val)) {
					FlatVector::SetNull(vec, row_count, true);
				} else {
					FlatVector::GetData<string_t>(vec)[row_count] = StringVector::AddString(vec, val);
				}
			}
		}

		row_count++;
	}

	output.SetCardinality(row_count);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPsamReader(ExtensionLoader &loader) {
	TableFunction read_psam("read_psam", {LogicalType::VARCHAR}, PsamScan, PsamBind, PsamInitGlobal, PsamInitLocal);
	read_psam.projection_pushdown = true;
	loader.RegisterFunction(read_psam);
}

} // namespace duckdb
