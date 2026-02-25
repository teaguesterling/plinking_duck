#include "psam_reader.hpp"
#include "duckdb.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/file_open_flags.hpp"

#include <algorithm>
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
// Line splitting utility
// ---------------------------------------------------------------------------

static vector<string> SplitTabLine(const string &line) {
	vector<string> fields;
	size_t start = 0;
	size_t pos = line.find('\t');
	while (pos != string::npos) {
		fields.push_back(line.substr(start, pos - start));
		start = pos + 1;
		pos = line.find('\t', start);
	}
	fields.push_back(line.substr(start));
	return fields;
}

// ---------------------------------------------------------------------------
// File reading via DuckDB VFS
// ---------------------------------------------------------------------------

//! Read an entire file via DuckDB's virtual file system and split into lines.
//! Strips \r from line endings. Returns an empty vector for empty files.
//! Uses VFS so that S3, HTTP, and custom filesystems work transparently.
static vector<string> ReadFileLines(ClientContext &context, const string &path) {
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	auto file_size = handle->GetFileSize();

	if (file_size == 0) {
		return {};
	}

	// Read entire file into a string buffer — .psam/.fam files are small
	// metadata files (typically KB to low MB), so this is safe and simple
	string content(file_size, '\0');
	handle->Read(const_cast<char *>(content.data()), file_size);

	// Split into lines, stripping \r and skipping trailing empty line
	vector<string> lines;
	size_t start = 0;
	for (size_t i = 0; i < content.size(); i++) {
		if (content[i] == '\n') {
			size_t end = i;
			if (end > start && content[end - 1] == '\r') {
				end--;
			}
			lines.push_back(content.substr(start, end - start));
			start = i + 1;
		}
	}
	// Handle last line without trailing newline
	if (start < content.size()) {
		size_t end = content.size();
		if (end > start && content[end - 1] == '\r') {
			end--;
		}
		if (end > start) {
			lines.push_back(content.substr(start, end - start));
		}
	}

	return lines;
}

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

PsamHeaderInfo ParsePsamHeader(ClientContext &context, const string &path) {
	auto lines = ReadFileLines(context, path);
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
		// Strip the leading # from the first field name
		first_line[0] = ' ';
		auto fields = SplitTabLine(first_line);
		// Trim leading whitespace from first field (was '#')
		if (!fields.empty()) {
			fields[0].erase(0, fields[0].find_first_not_of(' '));
		}

		if (fields.empty()) {
			throw IOException("read_psam: file '%s' has an empty header", path);
		}

		// Detect #FID vs #IID
		if (fields[0] == "FID") {
			info.format = PsamFormat::PSAM_FID;
		} else if (fields[0] == "IID") {
			info.format = PsamFormat::PSAM_IID;
		} else {
			throw IOException("read_psam: file '%s' header must start with #FID or #IID, got '#%s'",
			                  path, fields[0]);
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

// ---------------------------------------------------------------------------
// LoadSampleInfo — reusable utility for read_pgen / read_pfile
// ---------------------------------------------------------------------------

SampleInfo LoadSampleInfo(ClientContext &context, const string &path) {
	auto lines = ReadFileLines(context, path);
	if (lines.empty()) {
		throw IOException("read_psam: file '%s' is empty", path);
	}

	// Parse header to get format and column layout
	auto header = ParsePsamHeader(context, path);

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

		auto fields = SplitTabLine(line);
		if (iid_idx >= fields.size()) {
			throw IOException("read_psam: file '%s' line %d has %d fields, expected at least %d",
			                  path, i + 1, fields.size(), iid_idx + 1);
		}

		info.iids.push_back(fields[iid_idx]);
		if (has_fid) {
			info.fids.push_back(fields[fid_idx]);
		}
		info.iid_to_idx[fields[iid_idx]] = info.iids.size() - 1;
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
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PsamGlobalState : public GlobalTableFunctionState {
	//! All data lines from the file (pre-read during init for simplicity;
	//! .psam files are at most ~500MB at extreme biobank scale and typically
	//! much smaller — this is a metadata file, not genotype data)
	vector<vector<string>> rows;
	//! Next row to hand out to a scan call
	idx_t next_row_idx = 0;
	//! Mutex for thread-safe row claiming
	mutex lock;
	//! Projected column IDs (file column index for each output column)
	vector<column_t> column_ids;
	//! For each output column, whether it's a PAT/MAT column
	vector<bool> is_parent_col;

	idx_t MaxThreads() const override {
		return 1; // .psam files are small; single-threaded scan is fine
	}
};

// ---------------------------------------------------------------------------
// Local state
// ---------------------------------------------------------------------------

struct PsamLocalState : public LocalTableFunctionState {
	// No per-thread state needed for single-threaded scan
};

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

static unique_ptr<GlobalTableFunctionState> PsamInitGlobal(ClientContext &context,
                                                           TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PsamBindData>();
	auto state = make_uniq<PsamGlobalState>();

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

	auto lines = ReadFileLines(context, bind_data.file_path);
	idx_t expected_cols = bind_data.column_names.size();

	// Data lines start at index 1 for .psam (skip header), 0 for .fam
	idx_t data_start = (bind_data.format != PsamFormat::FAM) ? 1 : 0;

	for (idx_t i = data_start; i < lines.size(); i++) {
		auto &line = lines[i];
		if (line.empty()) {
			continue;
		}

		auto fields = SplitTabLine(line);

		if (fields.size() != expected_cols) {
			throw IOException(
			    "read_psam: file '%s' line %d has %d fields, expected %d",
			    bind_data.file_path, i + 1, fields.size(), expected_cols);
		}

		state->rows.push_back(std::move(fields));
	}

	return std::move(state);
}

static unique_ptr<LocalTableFunctionState> PsamInitLocal(ExecutionContext &context,
                                                         TableFunctionInitInput &input,
                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<PsamLocalState>();
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void PsamScan(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &bind_data = input.bind_data->Cast<PsamBindData>();
	auto &global_state = input.global_state->Cast<PsamGlobalState>();

	// Claim a batch of rows
	idx_t start_idx;
	idx_t batch_size;
	{
		lock_guard<mutex> guard(global_state.lock);
		start_idx = global_state.next_row_idx;
		batch_size = MinValue<idx_t>(STANDARD_VECTOR_SIZE,
		                             global_state.rows.size() - start_idx);
		global_state.next_row_idx += batch_size;
	}

	if (batch_size == 0) {
		output.SetCardinality(0);
		return;
	}

	auto &column_ids = global_state.column_ids;
	auto &is_parent_col = global_state.is_parent_col;

	for (idx_t row = 0; row < batch_size; row++) {
		auto &fields = global_state.rows[start_idx + row];

		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];

			// Handle ROW_ID pseudo-column
			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];
			const string &val = fields[file_col];

			if (file_col == bind_data.sex_col_idx) {
				// SEX column: integer with special missing handling
				// 1=male, 2=female; 0/NA/.  → NULL
				if (IsMissingValue(val)) {
					FlatVector::SetNull(vec, row, true);
				} else {
					try {
						int32_t sex_val = std::stoi(val);
						if (sex_val == 0) {
							FlatVector::SetNull(vec, row, true);
						} else {
							FlatVector::GetData<int32_t>(vec)[row] = sex_val;
						}
					} catch (...) {
						FlatVector::SetNull(vec, row, true);
					}
				}
			} else if (is_parent_col[out_col]) {
				// PAT/MAT columns: "0" means unknown parent → NULL
				if (IsParentMissing(val)) {
					FlatVector::SetNull(vec, row, true);
				} else {
					FlatVector::GetData<string_t>(vec)[row] = StringVector::AddString(vec, val);
				}
			} else {
				// VARCHAR columns: general missing-value handling
				if (IsMissingValue(val)) {
					FlatVector::SetNull(vec, row, true);
				} else {
					FlatVector::GetData<string_t>(vec)[row] = StringVector::AddString(vec, val);
				}
			}
		}
	}

	output.SetCardinality(batch_size);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPsamReader(ExtensionLoader &loader) {
	TableFunction read_psam("read_psam", {LogicalType::VARCHAR}, PsamScan, PsamBind, PsamInitGlobal,
	                        PsamInitLocal);
	read_psam.projection_pushdown = true;
	loader.RegisterFunction(read_psam);
}

} // namespace duckdb
