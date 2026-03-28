#include "plink_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

// ---------------------------------------------------------------------------
// Genotype output mode
// ---------------------------------------------------------------------------

GenotypeMode ResolveGenotypeMode(const string &mode_str, uint32_t sample_ct, const string &func_name) {
	auto mode = StringUtil::Lower(mode_str);
	if (mode == "auto") {
		return sample_ct <= ArrayType::MAX_ARRAY_SIZE ? GenotypeMode::ARRAY : GenotypeMode::LIST;
	} else if (mode == "array") {
		if (sample_ct > ArrayType::MAX_ARRAY_SIZE) {
			throw InvalidInputException("%s: genotypes := 'array' requires sample count (%u) <= %u. "
			                            "Use genotypes := 'list' or genotypes := 'auto' for large cohorts.",
			                            func_name, sample_ct, static_cast<uint32_t>(ArrayType::MAX_ARRAY_SIZE));
		}
		return GenotypeMode::ARRAY;
	} else if (mode == "list") {
		return GenotypeMode::LIST;
	} else if (mode == "columns") {
		return GenotypeMode::COLUMNS;
	} else if (mode == "struct") {
		return GenotypeMode::STRUCT;
	} else if (mode == "counts") {
		return GenotypeMode::COUNTS;
	} else if (mode == "stats") {
		return GenotypeMode::STATS;
	} else {
		throw InvalidInputException(
		    "%s: invalid genotypes value '%s' (expected 'auto', 'array', 'list', 'columns', 'struct', 'counts', or "
		    "'stats')",
		    func_name, mode_str);
	}
}

LogicalType MakeGenotypeCountsType() {
	child_list_t<LogicalType> children;
	children.push_back({"hom_ref", LogicalType::UINTEGER});
	children.push_back({"het", LogicalType::UINTEGER});
	children.push_back({"hom_alt", LogicalType::UINTEGER});
	children.push_back({"missing", LogicalType::UINTEGER});
	return LogicalType::STRUCT(std::move(children));
}

LogicalType MakeGenotypeStatsType() {
	child_list_t<LogicalType> children;
	children.push_back({"hom_ref", LogicalType::UINTEGER});
	children.push_back({"het", LogicalType::UINTEGER});
	children.push_back({"hom_alt", LogicalType::UINTEGER});
	children.push_back({"missing", LogicalType::UINTEGER});
	children.push_back({"n", LogicalType::UINTEGER});
	children.push_back({"af", LogicalType::DOUBLE});
	children.push_back({"maf", LogicalType::DOUBLE});
	children.push_back({"missing_rate", LogicalType::DOUBLE});
	children.push_back({"carrier_count", LogicalType::UINTEGER});
	children.push_back({"het_rate", LogicalType::DOUBLE});
	return LogicalType::STRUCT(std::move(children));
}

// ---------------------------------------------------------------------------
// Orient mode
// ---------------------------------------------------------------------------

OrientMode ResolveOrientMode(const string &orient_str, const string &func_name) {
	if (orient_str.empty()) {
		return OrientMode::VARIANT;
	}

	auto mode = StringUtil::Lower(orient_str);
	if (mode == "variant") {
		return OrientMode::VARIANT;
	} else if (mode == "genotype") {
		return OrientMode::GENOTYPE;
	} else if (mode == "sample") {
		return OrientMode::SAMPLE;
	} else {
		throw InvalidInputException("%s: invalid orient value '%s' (expected 'variant', 'genotype', or 'sample')",
		                            func_name, orient_str);
	}
}

// ---------------------------------------------------------------------------
// Offset-indexed variant metadata
// ---------------------------------------------------------------------------

size_t VariantMetadataIndex::LineEnd(idx_t vidx) const {
	size_t end;
	if (vidx + 1 < variant_ct) {
		end = static_cast<size_t>(line_offsets[vidx + 1]);
	} else {
		end = file_content.size();
	}
	// Strip trailing newline characters
	while (end > static_cast<size_t>(line_offsets[vidx]) &&
	       (file_content[end - 1] == '\n' || file_content[end - 1] == '\r')) {
		end--;
	}
	return end;
}

string VariantMetadataIndex::GetField(idx_t vidx, idx_t field_idx) const {
	auto start = static_cast<size_t>(line_offsets[vidx]);
	auto line_end = LineEnd(vidx);

	if (is_bim) {
		// Whitespace-delimited: skip leading whitespace, then find fields
		idx_t current_field = 0;
		size_t i = start;
		while (i < line_end) {
			// Skip whitespace
			while (i < line_end && (file_content[i] == ' ' || file_content[i] == '\t')) {
				i++;
			}
			if (i >= line_end) {
				break;
			}
			size_t field_start = i;
			// Find end of field
			while (i < line_end && file_content[i] != ' ' && file_content[i] != '\t') {
				i++;
			}
			if (current_field == field_idx) {
				return string(file_content.data() + field_start, i - field_start);
			}
			current_field++;
		}
	} else {
		// Tab-delimited
		idx_t current_field = 0;
		size_t field_start = start;
		for (size_t i = start; i < line_end; i++) {
			if (file_content[i] == '\t') {
				if (current_field == field_idx) {
					return string(file_content.data() + field_start, i - field_start);
				}
				current_field++;
				field_start = i + 1;
			}
		}
		// Last field on the line
		if (current_field == field_idx) {
			return string(file_content.data() + field_start, line_end - field_start);
		}
	}

	throw InternalException("VariantMetadataIndex::GetField: field index %llu out of range for variant %llu",
	                        static_cast<unsigned long long>(field_idx), static_cast<unsigned long long>(vidx));
}

string VariantMetadataIndex::GetChrom(idx_t vidx) const {
	return GetField(vidx, chrom_idx);
}

int32_t VariantMetadataIndex::GetPos(idx_t vidx) const {
	auto field = GetField(vidx, pos_idx);
	char *end;
	errno = 0;
	long val = std::strtol(field.c_str(), &end, 10);
	if (end == field.c_str() || *end != '\0' || errno != 0) {
		throw InternalException("VariantMetadataIndex::GetPos: invalid POS value '%s' for variant %llu", field.c_str(),
		                        static_cast<unsigned long long>(vidx));
	}
	return static_cast<int32_t>(val);
}

string VariantMetadataIndex::GetId(idx_t vidx) const {
	auto field = GetField(vidx, id_idx);
	if (field == ".") {
		return "";
	}
	return field;
}

string VariantMetadataIndex::GetRef(idx_t vidx) const {
	return GetField(vidx, ref_idx);
}

string VariantMetadataIndex::GetAlt(idx_t vidx) const {
	auto field = GetField(vidx, alt_idx);
	if (field == ".") {
		return "";
	}
	return field;
}

VariantMetadataIndex LoadVariantMetadataIndex(ClientContext &context, const string &path, const string &func_name) {
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	auto file_size = handle->GetFileSize();

	if (file_size == 0) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is empty", func_name, path);
	}

	VariantMetadataIndex idx;
	idx.file_content.resize(file_size);
	handle->Read(const_cast<char *>(idx.file_content.data()), file_size);

	// Parse header from the in-memory buffer (same logic as LoadVariantMetadata)
	size_t pos = 0;

	// Skip ## comment/meta lines
	while (pos < file_size) {
		if (file_size - pos >= 2 && idx.file_content[pos] == '#' && idx.file_content[pos + 1] == '#') {
			// Skip to next line
			while (pos < file_size && idx.file_content[pos] != '\n') {
				pos++;
			}
			if (pos < file_size) {
				pos++; // skip newline
			}
			continue;
		}
		if (pos < file_size && idx.file_content[pos] == '\n') {
			pos++;
			continue;
		}
		break;
	}

	if (pos >= file_size) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' contains no header or data", func_name, path);
	}

	// Extract the header/first-data line to determine format
	size_t header_start = pos;
	size_t header_end = pos;
	while (header_end < file_size && idx.file_content[header_end] != '\n') {
		header_end++;
	}
	// Strip trailing \r
	size_t header_content_end = header_end;
	if (header_content_end > header_start && idx.file_content[header_content_end - 1] == '\r') {
		header_content_end--;
	}
	string header_line(idx.file_content.data() + header_start, header_content_end - header_start);

	vector<string> column_names;

	if (header_line.size() >= 6 && header_line.substr(0, 6) == "#CHROM") {
		// .pvar format
		idx.is_bim = false;
		auto fields = SplitTabLine(header_line.substr(1)); // strip leading '#'
		column_names = std::move(fields);
		// Skip past header line
		pos = header_end;
		if (pos < file_size) {
			pos++; // skip newline
		}
	} else {
		// Legacy .bim format: no header to skip, first line is data
		idx.is_bim = true;
		// .bim physical columns: CHROM(0) ID(1) CM(2) POS(3) ALT(4) REF(5)
		// We store physical field indices directly
		column_names = {"CHROM", "ID", "CM", "POS", "ALT", "REF"};
	}

	// Find physical column indices
	idx.chrom_idx = DConstants::INVALID_INDEX;
	idx.pos_idx = DConstants::INVALID_INDEX;
	idx.id_idx = DConstants::INVALID_INDEX;
	idx.ref_idx = DConstants::INVALID_INDEX;
	idx.alt_idx = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < column_names.size(); i++) {
		const auto &name = column_names[i];
		if (name == "CHROM") {
			idx.chrom_idx = i;
		} else if (name == "POS") {
			idx.pos_idx = i;
		} else if (name == "ID") {
			idx.id_idx = i;
		} else if (name == "REF") {
			idx.ref_idx = i;
		} else if (name == "ALT") {
			idx.alt_idx = i;
		}
	}

	if (idx.chrom_idx == DConstants::INVALID_INDEX || idx.pos_idx == DConstants::INVALID_INDEX ||
	    idx.id_idx == DConstants::INVALID_INDEX || idx.ref_idx == DConstants::INVALID_INDEX ||
	    idx.alt_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is missing required columns "
		                            "(need CHROM, POS, ID, REF, ALT)",
		                            func_name, path);
	}

	// Build line offset index: one pass over the buffer
	// Estimate capacity to avoid reallocation
	if (file_size > 0) {
		idx.line_offsets.reserve(file_size / 30); // rough estimate: ~30 bytes per line
	}

	while (pos < file_size) {
		// Skip empty lines
		if (idx.file_content[pos] == '\n') {
			pos++;
			continue;
		}
		if (idx.file_content[pos] == '\r' && pos + 1 < file_size && idx.file_content[pos + 1] == '\n') {
			pos += 2;
			continue;
		}

		idx.line_offsets.push_back(static_cast<uint64_t>(pos));

		// Scan to next newline
		while (pos < file_size && idx.file_content[pos] != '\n') {
			pos++;
		}
		if (pos < file_size) {
			pos++; // skip newline
		}
	}

	idx.variant_ct = idx.line_offsets.size();
	return idx;
}

// ---------------------------------------------------------------------------
// File utilities
// ---------------------------------------------------------------------------

vector<string> ReadFileLines(ClientContext &context, const string &path) {
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	auto file_size = handle->GetFileSize();

	if (file_size == 0) {
		return {};
	}

	string content(file_size, '\0');
	handle->Read(const_cast<char *>(content.data()), file_size);

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

vector<string> SplitTabLine(const string &line) {
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

vector<string> SplitWhitespaceLine(const string &line) {
	vector<string> fields;
	size_t i = 0;
	while (i < line.size()) {
		while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) {
			i++;
		}
		if (i >= line.size()) {
			break;
		}
		size_t start = i;
		while (i < line.size() && line[i] != ' ' && line[i] != '\t') {
			i++;
		}
		fields.push_back(line.substr(start, i - start));
	}
	return fields;
}

// ---------------------------------------------------------------------------
// Companion file discovery
// ---------------------------------------------------------------------------

string ReplaceExtension(const string &path, const string &new_ext) {
	auto dot = path.rfind('.');
	if (dot == string::npos) {
		return path + new_ext;
	}
	return path.substr(0, dot) + new_ext;
}

string FindCompanionFile(FileSystem &fs, const string &pgen_path, const vector<string> &extensions) {
	for (auto &ext : extensions) {
		auto candidate = ReplaceExtension(pgen_path, ext);
		if (fs.FileExists(candidate)) {
			return candidate;
		}
	}
	return "";
}

bool IsParquetFile(const string &path) {
	return StringUtil::EndsWith(StringUtil::Lower(path), ".parquet");
}

bool IsNativePlinkFormat(const string &path) {
	auto lower = StringUtil::Lower(path);
	return StringUtil::EndsWith(lower, ".pvar") || StringUtil::EndsWith(lower, ".bim") ||
	       StringUtil::EndsWith(lower, ".psam") || StringUtil::EndsWith(lower, ".fam") ||
	       StringUtil::EndsWith(lower, ".pvar.zst") || StringUtil::EndsWith(lower, ".psam.zst");
}

bool GetUseParquetCompanions(ClientContext &context) {
	Value val;
	if (context.TryGetCurrentSetting("plinking_use_parquet_companions", val)) {
		return val.GetValue<bool>();
	}
	return true; // default
}

string FindCompanionFileWithParquet(ClientContext &context, FileSystem &fs, const string &pgen_path,
                                    const vector<string> &extensions) {
	// If parquet companions enabled, check for .parquet versions first
	if (GetUseParquetCompanions(context)) {
		for (auto &ext : extensions) {
			auto parquet_candidate = ReplaceExtension(pgen_path, ext + ".parquet");
			if (fs.FileExists(parquet_candidate)) {
				return parquet_candidate;
			}
		}
	}
	// Fall back to native text formats
	return FindCompanionFile(fs, pgen_path, extensions);
}

// ---------------------------------------------------------------------------
// Parquet companion loading
// ---------------------------------------------------------------------------

VariantMetadataIndex LoadVariantMetadataFromParquet(ClientContext &context, const string &path,
                                                    const string &func_name) {
	// Use a separate connection to read the parquet file (avoids bind reentrancy)
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto result = conn.TableFunction("parquet_scan", {Value(path)})->Execute();
	if (result->HasError()) {
		throw IOException("%s: failed to read parquet companion '%s': %s", func_name, path, result->GetError());
	}

	// Map columns by name (case-insensitive)
	auto &col_names = result->names;
	idx_t chrom_col = DConstants::INVALID_INDEX;
	idx_t pos_col = DConstants::INVALID_INDEX;
	idx_t id_col = DConstants::INVALID_INDEX;
	idx_t ref_col = DConstants::INVALID_INDEX;
	idx_t alt_col = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < col_names.size(); i++) {
		auto lower = StringUtil::Lower(col_names[i]);
		if (lower == "chrom" || lower == "#chrom") {
			chrom_col = i;
		} else if (lower == "pos") {
			pos_col = i;
		} else if (lower == "id") {
			id_col = i;
		} else if (lower == "ref") {
			ref_col = i;
		} else if (lower == "alt") {
			alt_col = i;
		}
	}

	if (chrom_col == DConstants::INVALID_INDEX || pos_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: parquet companion '%s' missing required columns "
		                            "(need CHROM and POS, found: %s)",
		                            func_name, path, StringUtil::Join(col_names, ", "));
	}

	// Synthesize a .pvar-format text buffer from the parquet data
	string buf;
	buf += "#CHROM\tPOS\tID\tREF\tALT\n";

	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		for (idx_t row = 0; row < chunk->size(); row++) {
			buf += chunk->GetValue(chrom_col, row).ToString();
			buf += '\t';
			buf += chunk->GetValue(pos_col, row).ToString();
			buf += '\t';
			buf += (id_col != DConstants::INVALID_INDEX) ? chunk->GetValue(id_col, row).ToString() : ".";
			buf += '\t';
			buf += (ref_col != DConstants::INVALID_INDEX) ? chunk->GetValue(ref_col, row).ToString() : ".";
			buf += '\t';
			buf += (alt_col != DConstants::INVALID_INDEX) ? chunk->GetValue(alt_col, row).ToString() : ".";
			buf += '\n';
		}
	}

	// Parse the synthetic buffer using the existing text parser
	VariantMetadataIndex idx;
	idx.file_content = std::move(buf);
	idx.is_bim = false;
	idx.chrom_idx = 0;
	idx.pos_idx = 1;
	idx.id_idx = 2;
	idx.ref_idx = 3;
	idx.alt_idx = 4;

	// Build line offsets — skip the header line (first line ending with \n)
	size_t pos = 0;
	// Skip header line
	while (pos < idx.file_content.size() && idx.file_content[pos] != '\n') {
		pos++;
	}
	if (pos < idx.file_content.size()) {
		pos++; // skip newline
	}

	// Index data lines
	while (pos < idx.file_content.size()) {
		if (idx.file_content[pos] == '\n') {
			pos++;
			continue;
		}
		idx.line_offsets.push_back(static_cast<uint64_t>(pos));
		while (pos < idx.file_content.size() && idx.file_content[pos] != '\n') {
			pos++;
		}
		if (pos < idx.file_content.size()) {
			pos++;
		}
	}

	idx.variant_ct = idx.line_offsets.size();
	return idx;
}

SampleInfo LoadSampleInfoFromParquet(ClientContext &context, const string &path) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto result = conn.TableFunction("parquet_scan", {Value(path)})->Execute();
	if (result->HasError()) {
		throw IOException("Failed to read parquet companion '%s': %s", path, result->GetError());
	}

	// Map columns by name (case-insensitive)
	auto &col_names = result->names;
	idx_t iid_col = DConstants::INVALID_INDEX;
	idx_t fid_col = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < col_names.size(); i++) {
		auto lower = StringUtil::Lower(col_names[i]);
		if (lower == "iid") {
			iid_col = i;
		} else if (lower == "fid") {
			fid_col = i;
		}
	}

	if (iid_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Parquet companion '%s' missing required IID column (found: %s)", path,
		                            StringUtil::Join(col_names, ", "));
	}

	SampleInfo info;
	bool has_fid = (fid_col != DConstants::INVALID_INDEX);

	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		for (idx_t row = 0; row < chunk->size(); row++) {
			auto iid = chunk->GetValue(iid_col, row).ToString();

			if (info.iid_to_idx.count(iid)) {
				throw IOException("Parquet companion '%s' has duplicate IID '%s'", path, iid);
			}

			info.iids.push_back(iid);
			if (has_fid) {
				info.fids.push_back(chunk->GetValue(fid_col, row).ToString());
			}
			info.iid_to_idx[iid] = info.iids.size() - 1;
		}
	}

	info.sample_ct = info.iids.size();
	return info;
}

// ---------------------------------------------------------------------------
// Generic source loading (CSV, tables, views, etc.)
// ---------------------------------------------------------------------------

//! Execute a query against a source string via a separate Connection.
//! The source may be a file path (CSV, parquet, etc.) or a table/view name.
//! DuckDB's replacement scan mechanism handles file extension dispatch.
static unique_ptr<MaterializedQueryResult> QuerySource(ClientContext &context, const string &source,
                                                       const string &func_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	// Escape single quotes for safety
	auto escaped = StringUtil::Replace(source, "'", "''");
	auto result = conn.Query("SELECT * FROM '" + escaped + "'");
	if (result->HasError()) {
		throw IOException("%s: failed to read source '%s': %s", func_name, source, result->GetError());
	}
	return result;
}

VariantMetadataIndex LoadVariantMetadataFromSource(ClientContext &context, const string &source,
                                                   const string &func_name) {
	auto result = QuerySource(context, source, func_name);

	// Map columns by name (case-insensitive) — same logic as parquet loader
	auto &col_names = result->names;
	idx_t chrom_col = DConstants::INVALID_INDEX;
	idx_t pos_col = DConstants::INVALID_INDEX;
	idx_t id_col = DConstants::INVALID_INDEX;
	idx_t ref_col = DConstants::INVALID_INDEX;
	idx_t alt_col = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < col_names.size(); i++) {
		auto lower = StringUtil::Lower(col_names[i]);
		if (lower == "chrom" || lower == "#chrom") {
			chrom_col = i;
		} else if (lower == "pos") {
			pos_col = i;
		} else if (lower == "id") {
			id_col = i;
		} else if (lower == "ref") {
			ref_col = i;
		} else if (lower == "alt") {
			alt_col = i;
		}
	}

	if (chrom_col == DConstants::INVALID_INDEX || pos_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: source '%s' missing required columns "
		                            "(need CHROM and POS, found: %s)",
		                            func_name, source, StringUtil::Join(col_names, ", "));
	}

	// Synthesize a .pvar-format text buffer from the query result
	string buf;
	buf += "#CHROM\tPOS\tID\tREF\tALT\n";

	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		for (idx_t row = 0; row < chunk->size(); row++) {
			buf += chunk->GetValue(chrom_col, row).ToString();
			buf += '\t';
			buf += chunk->GetValue(pos_col, row).ToString();
			buf += '\t';
			buf += (id_col != DConstants::INVALID_INDEX) ? chunk->GetValue(id_col, row).ToString() : ".";
			buf += '\t';
			buf += (ref_col != DConstants::INVALID_INDEX) ? chunk->GetValue(ref_col, row).ToString() : ".";
			buf += '\t';
			buf += (alt_col != DConstants::INVALID_INDEX) ? chunk->GetValue(alt_col, row).ToString() : ".";
			buf += '\n';
		}
	}

	// Parse the synthetic buffer using the existing text parser structure
	VariantMetadataIndex idx;
	idx.file_content = std::move(buf);
	idx.is_bim = false;
	idx.chrom_idx = 0;
	idx.pos_idx = 1;
	idx.id_idx = 2;
	idx.ref_idx = 3;
	idx.alt_idx = 4;

	// Build line offsets — skip the header line (first line ending with \n)
	size_t pos = 0;
	while (pos < idx.file_content.size() && idx.file_content[pos] != '\n') {
		pos++;
	}
	if (pos < idx.file_content.size()) {
		pos++; // skip newline
	}

	// Index data lines
	while (pos < idx.file_content.size()) {
		if (idx.file_content[pos] == '\n') {
			pos++;
			continue;
		}
		idx.line_offsets.push_back(static_cast<uint64_t>(pos));
		while (pos < idx.file_content.size() && idx.file_content[pos] != '\n') {
			pos++;
		}
		if (pos < idx.file_content.size()) {
			pos++;
		}
	}

	idx.variant_ct = idx.line_offsets.size();
	return idx;
}

SampleInfo LoadSampleInfoFromSource(ClientContext &context, const string &source) {
	auto result = QuerySource(context, source, "LoadSampleInfoFromSource");

	// Map columns by name (case-insensitive) — same logic as parquet loader
	auto &col_names = result->names;
	idx_t iid_col = DConstants::INVALID_INDEX;
	idx_t fid_col = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < col_names.size(); i++) {
		auto lower = StringUtil::Lower(col_names[i]);
		if (lower == "iid") {
			iid_col = i;
		} else if (lower == "fid") {
			fid_col = i;
		}
	}

	if (iid_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Source '%s' missing required IID column (found: %s)", source,
		                            StringUtil::Join(col_names, ", "));
	}

	SampleInfo info;
	bool has_fid = (fid_col != DConstants::INVALID_INDEX);

	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		for (idx_t row = 0; row < chunk->size(); row++) {
			auto iid = chunk->GetValue(iid_col, row).ToString();

			if (info.iid_to_idx.count(iid)) {
				throw IOException("Source '%s' has duplicate IID '%s'", source, iid);
			}

			info.iids.push_back(iid);
			if (has_fid) {
				info.fids.push_back(chunk->GetValue(fid_col, row).ToString());
			}
			info.iid_to_idx[iid] = info.iids.size() - 1;
		}
	}

	info.sample_ct = info.iids.size();
	return info;
}

// ---------------------------------------------------------------------------
// Unified dispatch functions
// ---------------------------------------------------------------------------

VariantMetadataIndex LoadVariantMetadata(ClientContext &context, const string &path, const string &func_name) {
	if (IsParquetFile(path)) {
		return LoadVariantMetadataFromParquet(context, path, func_name);
	}
	if (IsNativePlinkFormat(path)) {
		return LoadVariantMetadataIndex(context, path, func_name);
	}
	// Arbitrary source (CSV, table, view, etc.)
	return LoadVariantMetadataFromSource(context, path, func_name);
}

SampleInfo LoadSampleMetadata(ClientContext &context, const string &path) {
	if (IsParquetFile(path)) {
		return LoadSampleInfoFromParquet(context, path);
	}
	if (IsNativePlinkFormat(path)) {
		return LoadSampleInfo(context, path);
	}
	// Arbitrary source (CSV, table, view, etc.)
	return LoadSampleInfoFromSource(context, path);
}

// ---------------------------------------------------------------------------
// Sample parameter resolution
// ---------------------------------------------------------------------------

vector<uint32_t> ResolveSampleIndices(const Value &samples_val, uint32_t raw_sample_ct, const SampleInfo *sample_info,
                                      const string &func_name) {
	auto &child_type = ListType::GetChildType(samples_val.type());
	auto &children = ListValue::GetChildren(samples_val);

	if (children.empty()) {
		throw InvalidInputException("%s: samples list must not be empty", func_name);
	}

	vector<uint32_t> indices;

	if (child_type.id() == LogicalTypeId::INTEGER || child_type.id() == LogicalTypeId::BIGINT) {
		for (auto &child : children) {
			int64_t idx = child.GetValue<int64_t>();
			if (idx < 0 || static_cast<uint32_t>(idx) >= raw_sample_ct) {
				throw InvalidInputException("%s: sample index %lld out of range (sample count: %u)", func_name,
				                            static_cast<long long>(idx), raw_sample_ct);
			}
			indices.push_back(static_cast<uint32_t>(idx));
		}
	} else if (child_type.id() == LogicalTypeId::VARCHAR) {
		if (!sample_info) {
			throw InvalidInputException("%s: samples parameter requires LIST(INTEGER) when no .psam "
			                            "is available (no sample IDs to match against)",
			                            func_name);
		}
		for (auto &child : children) {
			auto iid = child.GetValue<string>();
			auto it = sample_info->iid_to_idx.find(iid);
			if (it == sample_info->iid_to_idx.end()) {
				throw InvalidInputException("%s: sample '%s' not found in .psam", func_name, iid);
			}
			indices.push_back(static_cast<uint32_t>(it->second));
		}
	} else {
		throw InvalidInputException("%s: samples parameter must be LIST(VARCHAR) or LIST(INTEGER)", func_name);
	}

	// Validate no duplicates
	std::unordered_set<uint32_t> seen;
	for (auto idx : indices) {
		if (!seen.insert(idx).second) {
			throw InvalidInputException("%s: duplicate sample index %u in samples list", func_name, idx);
		}
	}

	return indices;
}

// ---------------------------------------------------------------------------
// Sample subsetting
// ---------------------------------------------------------------------------

SampleSubset BuildSampleSubset(uint32_t raw_sample_ct, const vector<uint32_t> &sample_indices) {
	SampleSubset result;
	result.raw_sample_ct = raw_sample_ct;
	result.subset_sample_ct = static_cast<uint32_t>(sample_indices.size());

	// Allocate at vector-aligned size for FillInterleavedMaskVec
	uintptr_t aligned_word_ct = plink2::BitCtToAlignedWordCt(raw_sample_ct);
	uintptr_t alloc_size = aligned_word_ct * sizeof(uintptr_t);

	// sample_include bitmask
	result.sample_include_buf.Allocate(alloc_size);
	std::memset(result.sample_include_buf.ptr, 0, alloc_size);
	auto *sample_include = result.sample_include_buf.As<uintptr_t>();
	for (auto idx : sample_indices) {
		plink2::SetBit(idx, sample_include);
	}

	// Interleaved vec for PgrGetCounts
	result.interleaved_vec_buf.Allocate(alloc_size);
	uint32_t base_vec_ct = static_cast<uint32_t>(aligned_word_ct / plink2::kWordsPerVec);
	plink2::FillInterleavedMaskVec(sample_include, base_vec_ct, result.interleaved_vec_buf.As<uintptr_t>());

	// Cumulative popcounts for PgrSampleSubsetIndex
	uintptr_t word_ct = plink2::DivUp(raw_sample_ct, static_cast<uint32_t>(plink2::kBitsPerWord));
	result.cumulative_popcounts_buf.Allocate(word_ct * sizeof(uint32_t));
	plink2::FillCumulativePopcounts(sample_include, word_ct, result.cumulative_popcounts_buf.As<uint32_t>());

	return result;
}

// ---------------------------------------------------------------------------
// Region filtering
// ---------------------------------------------------------------------------

VariantRange ParseRegion(const string &region_str, const VariantMetadataIndex &variants, const string &func_name) {
	// Parse "chr:start-end"
	auto colon = region_str.find(':');
	if (colon == string::npos || colon == 0) {
		throw InvalidInputException("%s: invalid region format '%s' (expected 'chr:start-end')", func_name, region_str);
	}

	string chrom = region_str.substr(0, colon);
	string range_part = region_str.substr(colon + 1);

	auto dash = range_part.find('-');
	if (dash == string::npos) {
		throw InvalidInputException("%s: invalid region format '%s' (expected 'chr:start-end')", func_name, region_str);
	}

	string start_str = range_part.substr(0, dash);
	string end_str = range_part.substr(dash + 1);

	char *parse_end;
	errno = 0;
	long start_pos = std::strtol(start_str.c_str(), &parse_end, 10);
	if (parse_end == start_str.c_str() || *parse_end != '\0' || errno != 0 || start_pos < 0) {
		throw InvalidInputException("%s: invalid region start position in '%s'", func_name, region_str);
	}

	errno = 0;
	long end_pos = std::strtol(end_str.c_str(), &parse_end, 10);
	if (parse_end == end_str.c_str() || *parse_end != '\0' || errno != 0 || end_pos < 0) {
		throw InvalidInputException("%s: invalid region end position in '%s'", func_name, region_str);
	}

	// Scan variant metadata using on-demand field parsing
	VariantRange range;
	range.has_filter = true;

	bool found_start = false;
	for (uint32_t i = 0; i < static_cast<uint32_t>(variants.variant_ct); i++) {
		auto v_chrom = variants.GetField(i, variants.chrom_idx);
		if (v_chrom != chrom) {
			if (found_start) {
				break; // Past the matching chromosome block
			}
			continue;
		}
		auto v_pos = variants.GetPos(i);
		if (v_pos >= static_cast<int32_t>(start_pos) && v_pos <= static_cast<int32_t>(end_pos)) {
			if (!found_start) {
				range.start_idx = i;
				found_start = true;
			}
			range.end_idx = i + 1;
		} else if (found_start && v_pos > static_cast<int32_t>(end_pos)) {
			break; // Past the matching position range within the chromosome
		}
	}

	return range;
}

// ---------------------------------------------------------------------------
// Count-based filtering (af_range, ac_range)
// ---------------------------------------------------------------------------

RangeFilter ParseRangeFilter(const Value &val, const string &param_name, double valid_min, double valid_max,
                             const string &func_name) {
	RangeFilter result;

	if (val.type().id() != LogicalTypeId::STRUCT) {
		throw InvalidInputException("%s: %s must be a STRUCT (e.g. {min: 0.0, max: 0.5})", func_name, param_name);
	}

	auto &child_types = StructType::GetChildTypes(val.type());
	auto &children = StructValue::GetChildren(val);

	if (children.empty()) {
		// Empty struct = no filter
		return result;
	}

	for (idx_t i = 0; i < child_types.size(); i++) {
		auto &field_name = child_types[i].first;
		auto &child_val = children[i];

		if (field_name != "min" && field_name != "max") {
			throw InvalidInputException("%s: %s has unknown field '%s' (expected 'min' and/or 'max')", func_name,
			                            param_name, field_name);
		}

		if (child_val.IsNull()) {
			continue; // NULL = unbounded
		}

		double v = child_val.GetValue<double>();
		if (v < valid_min || v > valid_max) {
			throw InvalidInputException("%s: %s.%s value %g is out of range [%g, %g]", func_name, param_name,
			                            field_name, v, valid_min, valid_max);
		}

		if (field_name == "min") {
			result.min = v;
		} else {
			result.max = v;
		}
	}

	if (result.min > result.max) {
		throw InvalidInputException("%s: %s min (%g) > max (%g)", func_name, param_name, result.min, result.max);
	}

	result.active = true;
	return result;
}

bool VariantPassesCountFilter(const CountFilter &filter, const STD_ARRAY_REF(uint32_t, 4) genocounts,
                              uint32_t effective_sample_ct) {
	uint32_t non_missing = genocounts[0] + genocounts[1] + genocounts[2];
	if (non_missing == 0) {
		return false;
	}

	uint32_t ac = genocounts[1] + 2 * genocounts[2];

	if (filter.ac_filter.active && !filter.ac_filter.Passes(static_cast<double>(ac))) {
		return false;
	}

	if (filter.af_filter.active) {
		double af = static_cast<double>(ac) / (2.0 * static_cast<double>(non_missing));
		if (!filter.af_filter.Passes(af)) {
			return false;
		}
	}

	return true;
}

// ---------------------------------------------------------------------------
// Genotype range filtering (genotype_range)
// ---------------------------------------------------------------------------

GenotypeRangeResult CheckGenotypeRange(const RangeFilter &filter, const STD_ARRAY_REF(uint32_t, 4) genocounts) {
	GenotypeRangeResult result;
	result.any_pass = false;
	result.all_pass = true;

	// genocounts: [hom_ref(0), het(1), hom_alt(2), missing(3)]
	// Genotype values map directly to array indices 0, 1, 2
	for (int g = 0; g <= 2; g++) {
		bool in_range = filter.Passes(static_cast<double>(g));
		if (in_range && genocounts[g] > 0) {
			result.any_pass = true;
		}
		if (!in_range && genocounts[g] > 0) {
			result.all_pass = false;
		}
	}

	return result;
}

PreDecompFilterResult CheckPreDecompFilters(const CountFilter &count_filter, const GenotypeRangeFilter &genotype_filter,
                                            const STD_ARRAY_REF(uint32_t, 4) genocounts, uint32_t sample_ct) {
	PreDecompFilterResult result;
	result.skip = false;
	result.all_pass = true;

	if (count_filter.HasFilter() && !VariantPassesCountFilter(count_filter, genocounts, sample_ct)) {
		result.skip = true;
		return result;
	}

	if (genotype_filter.active) {
		auto gr = CheckGenotypeRange(genotype_filter.range, genocounts);
		if (!gr.any_pass) {
			result.skip = true;
			return result;
		}
		result.all_pass = gr.all_pass;
	}

	return result;
}

// ---------------------------------------------------------------------------
// Genotype normalization for PCA
// ---------------------------------------------------------------------------

VariantNorm ComputeVariantNorm(double alt_freq) {
	VariantNorm norm;
	if (alt_freq <= 0.0 || alt_freq >= 1.0) {
		norm.center = 0.0;
		norm.inv_stdev = 0.0;
		norm.skip = true;
		return norm;
	}
	norm.center = 2.0 * alt_freq;
	norm.inv_stdev = 1.0 / std::sqrt(2.0 * alt_freq * (1.0 - alt_freq));
	norm.skip = false;
	return norm;
}

void NormalizeGenotypes(const int8_t *genotypes, uint32_t sample_ct, const VariantNorm &norm, double *output) {
	for (uint32_t s = 0; s < sample_ct; s++) {
		if (genotypes[s] == -9) {
			output[s] = 0.0; // mean imputation (centered mean = 0)
		} else {
			output[s] = (static_cast<double>(genotypes[s]) - norm.center) * norm.inv_stdev;
		}
	}
}

// ---------------------------------------------------------------------------
// Phased genotype unpacking
// ---------------------------------------------------------------------------

void UnpackPhasedGenotypes(const int8_t *genotype_bytes, const uintptr_t *phasepresent, const uintptr_t *phaseinfo,
                           uint32_t sample_ct, int8_t *output_pairs) {
	for (uint32_t s = 0; s < sample_ct; s++) {
		int8_t geno = genotype_bytes[s];
		idx_t out_idx = static_cast<idx_t>(s) * 2;
		switch (geno) {
		case -9: // missing
			output_pairs[out_idx] = -9;
			output_pairs[out_idx + 1] = -9;
			break;
		case 0: // hom ref
			output_pairs[out_idx] = 0;
			output_pairs[out_idx + 1] = 0;
			break;
		case 2: // hom alt
			output_pairs[out_idx] = 1;
			output_pairs[out_idx + 1] = 1;
			break;
		case 1: // het — check phase
			if (plink2::IsSet(phasepresent, s) && plink2::IsSet(phaseinfo, s)) {
				// ALT|REF
				output_pairs[out_idx] = 1;
				output_pairs[out_idx + 1] = 0;
			} else {
				// REF|ALT (phased default) or unphased convention
				output_pairs[out_idx] = 0;
				output_pairs[out_idx + 1] = 1;
			}
			break;
		default: // shouldn't happen
			output_pairs[out_idx] = -9;
			output_pairs[out_idx + 1] = -9;
			break;
		}
	}
}

// ---------------------------------------------------------------------------
// Unified variants parameter resolution
// ---------------------------------------------------------------------------

//! Validate a 0-based variant index is within range.
static void ValidateVariantIndex(int64_t idx, uint32_t raw_variant_ct, const string &func_name) {
	if (idx < 0 || static_cast<uint32_t>(idx) >= raw_variant_ct) {
		throw InvalidInputException("%s: variant index %lld out of range (variant count: %u)", func_name,
		                            static_cast<long long>(idx), raw_variant_ct);
	}
}

//! Build a map from variant ID to 0-based index.
static unordered_map<string, uint32_t> BuildVariantIdIndex(const VariantMetadataIndex &variants) {
	unordered_map<string, uint32_t> id_to_idx;
	for (idx_t i = 0; i < variants.variant_ct; i++) {
		auto id = variants.GetId(i);
		if (!id.empty()) {
			id_to_idx[id] = static_cast<uint32_t>(i);
		}
	}
	return id_to_idx;
}

//! Resolve a single CPRA string (chrom:pos or chrom:pos:ref:alt) to a variant index.
static uint32_t ResolveCpraString(const string &cpra, const VariantMetadataIndex &variants, const string &func_name) {
	// Split on ':'
	vector<string> parts;
	size_t start = 0;
	for (size_t i = 0; i <= cpra.size(); i++) {
		if (i == cpra.size() || cpra[i] == ':') {
			parts.push_back(cpra.substr(start, i - start));
			start = i + 1;
		}
	}

	if (parts.size() != 2 && parts.size() != 4) {
		throw InvalidInputException("%s: invalid CPRA format '%s' (expected CHROM:POS or CHROM:POS:REF:ALT)", func_name,
		                            cpra);
	}

	auto &chrom = parts[0];
	char *end_ptr;
	errno = 0;
	long pos = std::strtol(parts[1].c_str(), &end_ptr, 10);
	if (*end_ptr != '\0' || errno != 0) {
		throw InvalidInputException("%s: invalid position in CPRA '%s'", func_name, cpra);
	}

	bool match_alleles = (parts.size() == 4);
	string ref_match, alt_match;
	if (match_alleles) {
		ref_match = parts[2];
		alt_match = parts[3];
	}

	for (idx_t i = 0; i < variants.variant_ct; i++) {
		if (variants.GetChrom(i) == chrom && variants.GetPos(i) == static_cast<int32_t>(pos)) {
			if (match_alleles) {
				if (variants.GetRef(i) == ref_match && variants.GetAlt(i) == alt_match) {
					return static_cast<uint32_t>(i);
				}
			} else {
				return static_cast<uint32_t>(i);
			}
		}
	}

	throw InvalidInputException("%s: variant '%s' not found", func_name, cpra);
}

//! Resolve a single variant string (rsid or CPRA) to a variant index.
static uint32_t ResolveVariantString(const string &id, const VariantMetadataIndex &variants,
                                     const unordered_map<string, uint32_t> &id_to_idx, const string &func_name) {
	// If it contains ':', treat as CPRA
	if (id.find(':') != string::npos) {
		return ResolveCpraString(id, variants, func_name);
	}

	// Otherwise, look up as rsid
	auto it = id_to_idx.find(id);
	if (it == id_to_idx.end()) {
		throw InvalidInputException("%s: variant '%s' not found", func_name, id);
	}
	return it->second;
}

//! Resolve a CPRA struct ({chrom, pos} or {chrom, pos, ref, alt}) to a variant index.
static uint32_t ResolveCpraStruct(const Value &val, const VariantMetadataIndex &variants, const string &func_name) {
	auto &child_types = StructType::GetChildTypes(val.type());
	auto &children = StructValue::GetChildren(val);

	string chrom;
	int32_t pos = 0;
	string ref_val, alt_val;
	bool has_ref = false, has_alt = false;

	for (idx_t i = 0; i < child_types.size(); i++) {
		auto &name = child_types[i].first;
		if (name == "chrom") {
			chrom = children[i].GetValue<string>();
		} else if (name == "pos") {
			pos = children[i].GetValue<int32_t>();
		} else if (name == "ref") {
			ref_val = children[i].GetValue<string>();
			has_ref = true;
		} else if (name == "alt") {
			alt_val = children[i].GetValue<string>();
			has_alt = true;
		}
	}

	bool match_alleles = has_ref && has_alt;

	for (idx_t i = 0; i < variants.variant_ct; i++) {
		if (variants.GetChrom(i) == chrom && variants.GetPos(i) == pos) {
			if (match_alleles) {
				if (variants.GetRef(i) == ref_val && variants.GetAlt(i) == alt_val) {
					return static_cast<uint32_t>(i);
				}
			} else {
				return static_cast<uint32_t>(i);
			}
		}
	}

	string desc = chrom + ":" + std::to_string(pos);
	if (match_alleles) {
		desc += ":" + ref_val + ":" + alt_val;
	}
	throw InvalidInputException("%s: variant '%s' not found", func_name, desc);
}

//! Resolve a range struct ({start, stop} with INTEGER or VARCHAR values) to variant indices.
static vector<uint32_t> ResolveRangeStruct(const Value &val, const VariantMetadataIndex &variants,
                                           uint32_t raw_variant_ct, const unordered_map<string, uint32_t> &id_to_idx,
                                           const string &func_name) {
	auto &child_types = StructType::GetChildTypes(val.type());
	auto &children = StructValue::GetChildren(val);

	idx_t start_field = DConstants::INVALID_INDEX;
	idx_t stop_field = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < child_types.size(); i++) {
		if (child_types[i].first == "start") {
			start_field = i;
		} else if (child_types[i].first == "stop") {
			stop_field = i;
		}
	}

	if (start_field == DConstants::INVALID_INDEX || stop_field == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: range struct must have 'start' and 'stop' fields", func_name);
	}

	auto &start_type = child_types[start_field].second;
	uint32_t start_idx, stop_idx;

	if (start_type.id() == LogicalTypeId::INTEGER || start_type.id() == LogicalTypeId::BIGINT) {
		int64_t sv = children[start_field].GetValue<int64_t>();
		int64_t ev = children[stop_field].GetValue<int64_t>();
		ValidateVariantIndex(sv, raw_variant_ct, func_name);
		ValidateVariantIndex(ev, raw_variant_ct, func_name);
		start_idx = static_cast<uint32_t>(sv);
		stop_idx = static_cast<uint32_t>(ev);
	} else if (start_type.id() == LogicalTypeId::VARCHAR) {
		auto start_str = children[start_field].GetValue<string>();
		auto stop_str = children[stop_field].GetValue<string>();
		start_idx = ResolveVariantString(start_str, variants, id_to_idx, func_name);
		stop_idx = ResolveVariantString(stop_str, variants, id_to_idx, func_name);
	} else {
		throw InvalidInputException("%s: range struct start/stop must be INTEGER or VARCHAR", func_name);
	}

	if (start_idx > stop_idx) {
		throw InvalidInputException("%s: variants range start (%u) is after stop (%u)", func_name, start_idx, stop_idx);
	}

	vector<uint32_t> indices;
	for (uint32_t i = start_idx; i <= stop_idx; i++) {
		indices.push_back(i);
	}
	return indices;
}

vector<uint32_t> ResolveVariantsParameter(const Value &val, const VariantMetadataIndex &variants,
                                          uint32_t raw_variant_ct, const string &func_name) {
	vector<uint32_t> indices;
	auto &type = val.type();

	// Lazily build ID index only when needed
	unordered_map<string, uint32_t> id_to_idx;
	auto ensure_id_index = [&]() {
		if (id_to_idx.empty()) {
			id_to_idx = BuildVariantIdIndex(variants);
		}
	};

	if (type.id() == LogicalTypeId::INTEGER || type.id() == LogicalTypeId::BIGINT) {
		// Single integer index
		int64_t idx = val.GetValue<int64_t>();
		ValidateVariantIndex(idx, raw_variant_ct, func_name);
		indices.push_back(static_cast<uint32_t>(idx));

	} else if (type.id() == LogicalTypeId::VARCHAR) {
		// Single rsid or CPRA string
		ensure_id_index();
		indices.push_back(ResolveVariantString(val.GetValue<string>(), variants, id_to_idx, func_name));

	} else if (type.id() == LogicalTypeId::STRUCT) {
		// Could be CPRA struct or range struct — disambiguate by field names
		auto &child_types = StructType::GetChildTypes(type);
		bool has_start = false, has_chrom = false;
		for (auto &ct : child_types) {
			if (ct.first == "start") {
				has_start = true;
			}
			if (ct.first == "chrom") {
				has_chrom = true;
			}
		}

		if (has_start && has_chrom) {
			throw InvalidInputException("%s: ambiguous variants struct — has both 'start' and 'chrom' fields. "
			                            "Use {start:, stop:} for a range or {chrom:, pos:} for a CPRA lookup.",
			                            func_name);
		} else if (has_start) {
			ensure_id_index();
			indices = ResolveRangeStruct(val, variants, raw_variant_ct, id_to_idx, func_name);
		} else if (has_chrom) {
			indices.push_back(ResolveCpraStruct(val, variants, func_name));
		} else {
			throw InvalidInputException(
			    "%s: variants struct must have either 'start'/'stop' (range) or 'chrom'/'pos' (CPRA) fields",
			    func_name);
		}

	} else if (type.id() == LogicalTypeId::LIST) {
		auto &child_type = ListType::GetChildType(type);
		auto &children = ListValue::GetChildren(val);

		if (children.empty()) {
			throw InvalidInputException("%s: variants list must not be empty", func_name);
		}

		if (child_type.id() == LogicalTypeId::INTEGER || child_type.id() == LogicalTypeId::BIGINT) {
			for (auto &child : children) {
				int64_t idx = child.GetValue<int64_t>();
				ValidateVariantIndex(idx, raw_variant_ct, func_name);
				indices.push_back(static_cast<uint32_t>(idx));
			}
		} else if (child_type.id() == LogicalTypeId::VARCHAR) {
			ensure_id_index();
			for (auto &child : children) {
				auto id = child.GetValue<string>();
				indices.push_back(ResolveVariantString(id, variants, id_to_idx, func_name));
			}
		} else if (child_type.id() == LogicalTypeId::STRUCT) {
			for (auto &child : children) {
				indices.push_back(ResolveCpraStruct(child, variants, func_name));
			}
		} else {
			throw InvalidInputException("%s: variants list elements must be INTEGER, VARCHAR, or STRUCT (got %s)",
			                            func_name, child_type.ToString());
		}

	} else {
		throw InvalidInputException("%s: variants parameter must be an integer, string, struct, or list (got %s)",
		                            func_name, type.ToString());
	}

	// Validate no duplicates
	{
		std::unordered_set<uint32_t> seen;
		for (auto idx : indices) {
			if (!seen.insert(idx).second) {
				throw InvalidInputException("%s: duplicate variant index %u in variants parameter", func_name, idx);
			}
		}
	}

	return indices;
}

// ---------------------------------------------------------------------------
// Max threads config helper
// ---------------------------------------------------------------------------

uint32_t GetPlinkingMaxThreads(ClientContext &context) {
	Value val;
	if (context.TryGetCurrentSetting("plinking_max_threads", val)) {
		auto v = val.GetValue<int64_t>();
		if (v > 0) {
			return static_cast<uint32_t>(v);
		}
	}
	return 0;
}

idx_t ApplyMaxThreadsCap(idx_t computed, uint32_t config_max_threads) {
	if (config_max_threads > 0) {
		return MinValue<idx_t>(computed, static_cast<idx_t>(config_max_threads));
	}
	return MinValue<idx_t>(computed, 16);
}

} // namespace duckdb
