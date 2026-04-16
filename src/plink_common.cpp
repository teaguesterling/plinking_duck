#include "plink_common.hpp"
#include "plink_profile.hpp"

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

// ---------------------------------------------------------------------------
// Columnar load helpers (shared between parquet and text paths)
// ---------------------------------------------------------------------------

//! Build chrom_offsets from the chroms column. Assumes (CHROM, POS)-sorted input.
static void BuildChromOffsets(VariantMetadataIndex &idx) {
	idx.chrom_offsets.clear();
	if (idx.chroms.empty()) {
		return;
	}
	idx_t run_start = 0;
	const string *current = &idx.chroms[0];
	for (idx_t i = 1; i < idx.chroms.size(); i++) {
		if (idx.chroms[i] != *current) {
			idx.chrom_offsets.emplace(*current, std::make_pair(run_start, i));
			current = &idx.chroms[i];
			run_start = i;
		}
	}
	idx.chrom_offsets.emplace(*current, std::make_pair(run_start, idx.chroms.size()));
}

//! Parse a single delimited field from [start, end) in buf (no allocation).
//! Advances *field_start to the start of the next field (after the delimiter),
//! sets *out_start and *out_len to the span of the current field within buf.
//! Returns false if no more fields are found. Delim = '\t' for pvar, ' '/'\t' for bim.
static bool NextField(const char *buf, size_t line_end, size_t *cursor, size_t *out_start, size_t *out_len,
                      bool whitespace) {
	size_t i = *cursor;
	if (whitespace) {
		while (i < line_end && (buf[i] == ' ' || buf[i] == '\t')) {
			i++;
		}
		if (i >= line_end) {
			return false;
		}
		size_t start = i;
		while (i < line_end && buf[i] != ' ' && buf[i] != '\t') {
			i++;
		}
		*out_start = start;
		*out_len = i - start;
		*cursor = i;
	} else {
		if (i > line_end) {
			return false;
		}
		size_t start = i;
		while (i < line_end && buf[i] != '\t') {
			i++;
		}
		*out_start = start;
		*out_len = i - start;
		*cursor = i < line_end ? i + 1 : i; // step past the tab
	}
	return true;
}

VariantMetadataIndex LoadVariantMetadataIndex(ClientContext &context, const string &path, const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataIndex(text)");
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	auto file_size = handle->GetFileSize();
	timer.Note("file_size=%llu", static_cast<unsigned long long>(file_size));

	if (file_size == 0) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is empty", func_name, path);
	}

	string file_content;
	file_content.resize(file_size);
	handle->Read(const_cast<char *>(file_content.data()), file_size);

	VariantMetadataIndex idx;

	// --- Header / column layout ---
	size_t pos = 0;
	while (pos < file_size) {
		if (file_size - pos >= 2 && file_content[pos] == '#' && file_content[pos + 1] == '#') {
			while (pos < file_size && file_content[pos] != '\n') {
				pos++;
			}
			if (pos < file_size) {
				pos++;
			}
			continue;
		}
		if (pos < file_size && file_content[pos] == '\n') {
			pos++;
			continue;
		}
		break;
	}

	if (pos >= file_size) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' contains no header or data", func_name, path);
	}

	size_t header_start = pos;
	size_t header_end = pos;
	while (header_end < file_size && file_content[header_end] != '\n') {
		header_end++;
	}
	size_t header_content_end = header_end;
	if (header_content_end > header_start && file_content[header_content_end - 1] == '\r') {
		header_content_end--;
	}
	string header_line(file_content.data() + header_start, header_content_end - header_start);

	idx_t chrom_field = DConstants::INVALID_INDEX;
	idx_t pos_field = DConstants::INVALID_INDEX;
	idx_t id_field = DConstants::INVALID_INDEX;
	idx_t ref_field = DConstants::INVALID_INDEX;
	idx_t alt_field = DConstants::INVALID_INDEX;

	if (header_line.size() >= 6 && header_line.substr(0, 6) == "#CHROM") {
		idx.is_bim = false;
		auto fields = SplitTabLine(header_line.substr(1));
		for (idx_t i = 0; i < fields.size(); i++) {
			if (fields[i] == "CHROM") {
				chrom_field = i;
			} else if (fields[i] == "POS") {
				pos_field = i;
			} else if (fields[i] == "ID") {
				id_field = i;
			} else if (fields[i] == "REF") {
				ref_field = i;
			} else if (fields[i] == "ALT") {
				alt_field = i;
			}
		}
		pos = header_end;
		if (pos < file_size) {
			pos++;
		}
	} else {
		idx.is_bim = true;
		// .bim: CHROM(0) ID(1) CM(2) POS(3) ALT(4) REF(5)
		chrom_field = 0;
		id_field = 1;
		pos_field = 3;
		alt_field = 4;
		ref_field = 5;
	}

	if (chrom_field == DConstants::INVALID_INDEX || pos_field == DConstants::INVALID_INDEX ||
	    id_field == DConstants::INVALID_INDEX || ref_field == DConstants::INVALID_INDEX ||
	    alt_field == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is missing required columns "
		                            "(need CHROM, POS, ID, REF, ALT)",
		                            func_name, path);
	}

	// --- Pre-size vectors (rough estimate) ---
	idx_t capacity = file_size / 30;
	idx.chroms.reserve(capacity);
	idx.positions.reserve(capacity);
	idx.ids.reserve(capacity);
	idx.refs.reserve(capacity);
	idx.alts.reserve(capacity);

	// --- Parse each data line directly into columnar vectors ---
	const idx_t max_field = std::max({chrom_field, pos_field, id_field, ref_field, alt_field});
	const char *buf = file_content.data();
	const bool whitespace = idx.is_bim;

	while (pos < file_size) {
		if (buf[pos] == '\n' || (buf[pos] == '\r' && pos + 1 < file_size && buf[pos + 1] == '\n')) {
			pos += (buf[pos] == '\r') ? 2 : 1;
			continue;
		}
		size_t line_end = pos;
		while (line_end < file_size && buf[line_end] != '\n') {
			line_end++;
		}
		size_t content_end = line_end;
		if (content_end > pos && buf[content_end - 1] == '\r') {
			content_end--;
		}

		size_t cursor = pos;
		size_t fstart = 0, flen = 0;
		string chrom, id, ref, alt;
		int32_t pos_val = 0;
		for (idx_t f = 0; f <= max_field; f++) {
			if (!NextField(buf, content_end, &cursor, &fstart, &flen, whitespace)) {
				break;
			}
			if (f == chrom_field) {
				chrom.assign(buf + fstart, flen);
			} else if (f == pos_field) {
				// strtol on a non-nul-terminated span: copy into small local
				char tmp[32];
				idx_t n = flen < sizeof(tmp) - 1 ? flen : sizeof(tmp) - 1;
				std::memcpy(tmp, buf + fstart, n);
				tmp[n] = '\0';
				char *end;
				errno = 0;
				long v = std::strtol(tmp, &end, 10);
				if (end == tmp || *end != '\0' || errno != 0) {
					throw InvalidInputException("%s: invalid POS value '%s' at line offset %llu", func_name, tmp,
					                            static_cast<unsigned long long>(pos));
				}
				pos_val = static_cast<int32_t>(v);
			} else if (f == id_field) {
				if (flen == 1 && buf[fstart] == '.') {
					// empty
				} else {
					id.assign(buf + fstart, flen);
				}
			} else if (f == ref_field) {
				ref.assign(buf + fstart, flen);
			} else if (f == alt_field) {
				if (flen == 1 && buf[fstart] == '.') {
					// empty
				} else {
					alt.assign(buf + fstart, flen);
				}
			}
		}

		idx.chroms.emplace_back(std::move(chrom));
		idx.positions.push_back(pos_val);
		idx.ids.emplace_back(std::move(id));
		idx.refs.emplace_back(std::move(ref));
		idx.alts.emplace_back(std::move(alt));

		pos = line_end < file_size ? line_end + 1 : line_end;
	}

	idx.variant_ct = idx.chroms.size();
	idx.has_ids = true;
	idx.has_alleles = true;
	BuildChromOffsets(idx);
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

//! Columnar ingestion helper. Scans a QueryResult and populates idx with
//! columnar-typed vectors. Maps columns by name (case-insensitive).
//! row_number_col: if != INVALID_INDEX, interpreted as a BIGINT file_row_number
//! column; the resulting index is sparse (idx.vidx_map populated). Otherwise
//! the index is dense (indexed by emit order).
static void IngestVariantResult(QueryResult &result, const string &source_label, const string &func_name,
                                VariantMetadataIndex &idx) {
	auto &col_names = result.names;
	auto &col_types = result.types;
	idx_t chrom_col = DConstants::INVALID_INDEX;
	idx_t pos_col = DConstants::INVALID_INDEX;
	idx_t id_col = DConstants::INVALID_INDEX;
	idx_t ref_col = DConstants::INVALID_INDEX;
	idx_t alt_col = DConstants::INVALID_INDEX;
	idx_t rn_col = DConstants::INVALID_INDEX;

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
		} else if (lower == "file_row_number") {
			rn_col = i;
		}
	}

	if (chrom_col == DConstants::INVALID_INDEX || pos_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: %s missing required columns (need CHROM and POS, found: %s)", func_name,
		                            source_label, StringUtil::Join(col_names, ", "));
	}

	auto pos_type_id = col_types[pos_col].id();
	if (pos_type_id != LogicalTypeId::INTEGER && pos_type_id != LogicalTypeId::BIGINT &&
	    pos_type_id != LogicalTypeId::UINTEGER && pos_type_id != LogicalTypeId::UBIGINT) {
		throw InvalidInputException("%s: %s POS column has unsupported type %s (expected INTEGER/BIGINT)", func_name,
		                            source_label, col_types[pos_col].ToString());
	}
	auto chrom_type_id = col_types[chrom_col].id();
	bool chrom_is_integer = (chrom_type_id == LogicalTypeId::INTEGER || chrom_type_id == LogicalTypeId::BIGINT ||
	                         chrom_type_id == LogicalTypeId::UINTEGER || chrom_type_id == LogicalTypeId::UBIGINT);

	unique_ptr<DataChunk> chunk;
	while ((chunk = result.Fetch()) != nullptr && chunk->size() > 0) {
		idx_t n = chunk->size();

		// CHROM — may be VARCHAR or numeric
		{
			auto &vec = chunk->data[chrom_col];
			UnifiedVectorFormat uvf;
			vec.ToUnifiedFormat(n, uvf);
			if (chrom_is_integer) {
				// Render numeric chrom as string for consistency with text .pvar
				for (idx_t i = 0; i < n; i++) {
					auto si = uvf.sel->get_index(i);
					if (!uvf.validity.RowIsValid(si)) {
						idx.chroms.emplace_back();
					} else {
						int64_t v;
						if (chrom_type_id == LogicalTypeId::INTEGER) {
							v = reinterpret_cast<const int32_t *>(uvf.data)[si];
						} else if (chrom_type_id == LogicalTypeId::BIGINT) {
							v = reinterpret_cast<const int64_t *>(uvf.data)[si];
						} else if (chrom_type_id == LogicalTypeId::UINTEGER) {
							v = reinterpret_cast<const uint32_t *>(uvf.data)[si];
						} else {
							v = static_cast<int64_t>(reinterpret_cast<const uint64_t *>(uvf.data)[si]);
						}
						idx.chroms.emplace_back(std::to_string(v));
					}
				}
			} else {
				auto data = reinterpret_cast<const string_t *>(uvf.data);
				for (idx_t i = 0; i < n; i++) {
					auto si = uvf.sel->get_index(i);
					if (!uvf.validity.RowIsValid(si)) {
						idx.chroms.emplace_back();
					} else {
						idx.chroms.emplace_back(data[si].GetData(), data[si].GetSize());
					}
				}
			}
		}

		// POS — widen to int32 (pgen uses int32)
		{
			auto &vec = chunk->data[pos_col];
			UnifiedVectorFormat uvf;
			vec.ToUnifiedFormat(n, uvf);
			for (idx_t i = 0; i < n; i++) {
				auto si = uvf.sel->get_index(i);
				int32_t v = 0;
				if (uvf.validity.RowIsValid(si)) {
					if (pos_type_id == LogicalTypeId::INTEGER) {
						v = reinterpret_cast<const int32_t *>(uvf.data)[si];
					} else if (pos_type_id == LogicalTypeId::BIGINT) {
						v = static_cast<int32_t>(reinterpret_cast<const int64_t *>(uvf.data)[si]);
					} else if (pos_type_id == LogicalTypeId::UINTEGER) {
						v = static_cast<int32_t>(reinterpret_cast<const uint32_t *>(uvf.data)[si]);
					} else {
						v = static_cast<int32_t>(reinterpret_cast<const uint64_t *>(uvf.data)[si]);
					}
				}
				idx.positions.push_back(v);
			}
		}

		// ID — VARCHAR, "." → ""
		{
			if (id_col == DConstants::INVALID_INDEX) {
				for (idx_t i = 0; i < n; i++) {
					idx.ids.emplace_back();
				}
			} else {
				auto &vec = chunk->data[id_col];
				UnifiedVectorFormat uvf;
				vec.ToUnifiedFormat(n, uvf);
				auto data = reinterpret_cast<const string_t *>(uvf.data);
				for (idx_t i = 0; i < n; i++) {
					auto si = uvf.sel->get_index(i);
					if (!uvf.validity.RowIsValid(si)) {
						idx.ids.emplace_back();
					} else {
						auto sz = data[si].GetSize();
						auto d = data[si].GetData();
						if (sz == 1 && d[0] == '.') {
							idx.ids.emplace_back();
						} else {
							idx.ids.emplace_back(d, sz);
						}
					}
				}
			}
		}

		// REF — VARCHAR (plain copy)
		{
			if (ref_col == DConstants::INVALID_INDEX) {
				for (idx_t i = 0; i < n; i++) {
					idx.refs.emplace_back();
				}
			} else {
				auto &vec = chunk->data[ref_col];
				UnifiedVectorFormat uvf;
				vec.ToUnifiedFormat(n, uvf);
				auto data = reinterpret_cast<const string_t *>(uvf.data);
				for (idx_t i = 0; i < n; i++) {
					auto si = uvf.sel->get_index(i);
					if (!uvf.validity.RowIsValid(si)) {
						idx.refs.emplace_back();
					} else {
						idx.refs.emplace_back(data[si].GetData(), data[si].GetSize());
					}
				}
			}
		}

		// ALT — VARCHAR, "." → ""
		{
			if (alt_col == DConstants::INVALID_INDEX) {
				for (idx_t i = 0; i < n; i++) {
					idx.alts.emplace_back();
				}
			} else {
				auto &vec = chunk->data[alt_col];
				UnifiedVectorFormat uvf;
				vec.ToUnifiedFormat(n, uvf);
				auto data = reinterpret_cast<const string_t *>(uvf.data);
				for (idx_t i = 0; i < n; i++) {
					auto si = uvf.sel->get_index(i);
					if (!uvf.validity.RowIsValid(si)) {
						idx.alts.emplace_back();
					} else {
						auto sz = data[si].GetSize();
						auto d = data[si].GetData();
						if (sz == 1 && d[0] == '.') {
							idx.alts.emplace_back();
						} else {
							idx.alts.emplace_back(d, sz);
						}
					}
				}
			}
		}

		// file_row_number → sparse vidx_map
		if (rn_col != DConstants::INVALID_INDEX) {
			auto &vec = chunk->data[rn_col];
			UnifiedVectorFormat uvf;
			vec.ToUnifiedFormat(n, uvf);
			auto data = reinterpret_cast<const int64_t *>(uvf.data);
			idx_t local_base = idx.chroms.size() - n;
			for (idx_t i = 0; i < n; i++) {
				auto si = uvf.sel->get_index(i);
				if (!uvf.validity.RowIsValid(si)) {
					throw IOException("%s: NULL file_row_number encountered (internal error)", func_name);
				}
				idx.vidx_map.emplace(static_cast<uint32_t>(data[si]), static_cast<uint32_t>(local_base + i));
			}
		}
	}
}

VariantMetadataIndex LoadVariantMetadataFromParquet(ClientContext &context, const string &path,
                                                    const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataFromParquet");
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto result = conn.TableFunction("parquet_scan", {Value(path)})->Execute();
	timer.Note("parquet_scan opened");
	if (result->HasError()) {
		throw IOException("%s: failed to read parquet companion '%s': %s", func_name, path, result->GetError());
	}

	VariantMetadataIndex idx;
	idx.is_bim = false;
	string label = "parquet companion '" + path + "'";
	IngestVariantResult(*result, label, func_name, idx);
	idx.variant_ct = idx.chroms.size();
	idx.has_ids = true;
	idx.has_alleles = true;
	timer.Note("ingested %llu variants", (unsigned long long)idx.variant_ct);
	BuildChromOffsets(idx);
	timer.Note("built chrom_offsets (%llu chroms)", (unsigned long long)idx.chrom_offsets.size());
	return idx;
}

//! Region-pushdown loader: queries the parquet file with a WHERE clause so
//! only matching rows are materialized. Returns a sparse index whose
//! vidx_map keys are the file row numbers (pgenlib-compatible vidx).
//! total_row_ct lets caller set variant_ct correctly for count validation.
VariantMetadataIndex LoadVariantMetadataFromParquetRegion(ClientContext &context, const string &path,
                                                          const string &chrom, int64_t pos_start, int64_t pos_end,
                                                          idx_t total_row_ct, const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataFromParquetRegion");
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto escaped_path = StringUtil::Replace(path, "'", "''");
	auto escaped_chrom = StringUtil::Replace(chrom, "'", "''");

	// file_row_number gives pgenlib-compatible global vidx; stats pruning + column
	// projection keep this O(region_size), not O(total).
	string sql = "SELECT \"#CHROM\" AS \"#CHROM\", POS, ID, REF, ALT, file_row_number "
	             "FROM read_parquet('" +
	             escaped_path +
	             "', file_row_number=true) "
	             "WHERE (CAST(\"#CHROM\" AS VARCHAR) = '" +
	             escaped_chrom + "') AND POS BETWEEN " + std::to_string(pos_start) + " AND " + std::to_string(pos_end) +
	             " ORDER BY file_row_number";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		// Retry without `#` prefix column alias if the file uses "CHROM" (alias still selects by position)
		auto alt_sql = "SELECT CHROM AS \"#CHROM\", POS, ID, REF, ALT, file_row_number "
		               "FROM read_parquet('" +
		               escaped_path +
		               "', file_row_number=true) "
		               "WHERE (CAST(CHROM AS VARCHAR) = '" +
		               escaped_chrom + "') AND POS BETWEEN " + std::to_string(pos_start) + " AND " +
		               std::to_string(pos_end) + " ORDER BY file_row_number";
		result = conn.Query(alt_sql);
		if (result->HasError()) {
			throw IOException("%s: region-pushdown query failed on '%s': %s", func_name, path, result->GetError());
		}
	}
	timer.Note("pushdown query returned");

	VariantMetadataIndex idx;
	idx.is_bim = false;
	string label = "parquet companion '" + path + "' (region pushdown)";
	IngestVariantResult(*result, label, func_name, idx);
	idx.variant_ct = total_row_ct; // total from parquet metadata (for count-mismatch validation)
	idx.has_ids = true;
	idx.has_alleles = true;
	// Single chrom in pushdown, single contiguous range of local indices
	if (!idx.chroms.empty()) {
		idx.chrom_offsets.emplace(idx.chroms.front(), std::make_pair(idx_t {0}, idx.chroms.size()));
	}
	timer.Note("ingested %llu region variants (total_ct=%llu)", (unsigned long long)idx.chroms.size(),
	           (unsigned long long)total_row_ct);
	return idx;
}

//! Get total row count from a parquet file via its footer metadata (no row scan).
idx_t GetParquetRowCount(ClientContext &context, const string &path) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto escaped = StringUtil::Replace(path, "'", "''");
	auto result = conn.Query("SELECT COUNT(*) FROM read_parquet('" + escaped + "')");
	if (result->HasError()) {
		throw IOException("Failed to get row count for '%s': %s", path, result->GetError());
	}
	auto chunk = result->Fetch();
	if (!chunk || chunk->size() == 0) {
		return 0;
	}
	return static_cast<idx_t>(chunk->GetValue(0, 0).GetValue<int64_t>());
}

//! Columnar ingest helper for psam query results.
static void IngestSampleResult(QueryResult &result, const string &source_label, SampleInfo &info, bool load_iids,
                               bool load_fids) {
	auto &col_names = result.names;
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
		throw InvalidInputException("%s missing required IID column (found: %s)", source_label,
		                            StringUtil::Join(col_names, ", "));
	}
	bool has_fid_col = (fid_col != DConstants::INVALID_INDEX);

	idx_t total_rows = 0;
	unique_ptr<DataChunk> chunk;
	while ((chunk = result.Fetch()) != nullptr && chunk->size() > 0) {
		idx_t n = chunk->size();
		total_rows += n;

		if (load_iids) {
			auto &vec = chunk->data[iid_col];
			UnifiedVectorFormat uvf;
			vec.ToUnifiedFormat(n, uvf);
			auto data = reinterpret_cast<const string_t *>(uvf.data);
			for (idx_t i = 0; i < n; i++) {
				auto si = uvf.sel->get_index(i);
				if (!uvf.validity.RowIsValid(si)) {
					info.iids.emplace_back();
				} else {
					info.iids.emplace_back(data[si].GetData(), data[si].GetSize());
				}
			}
		}
		if (load_fids && has_fid_col) {
			auto &vec = chunk->data[fid_col];
			UnifiedVectorFormat uvf;
			vec.ToUnifiedFormat(n, uvf);
			auto data = reinterpret_cast<const string_t *>(uvf.data);
			for (idx_t i = 0; i < n; i++) {
				auto si = uvf.sel->get_index(i);
				if (!uvf.validity.RowIsValid(si)) {
					info.fids.emplace_back();
				} else {
					info.fids.emplace_back(data[si].GetData(), data[si].GetSize());
				}
			}
		}
	}
	if (!load_iids) {
		info.sample_ct = total_rows;
	}
}

SampleInfo LoadSampleInfoFromParquet(ClientContext &context, const string &path) {
	BindPhaseTimer timer("LoadSampleInfoFromParquet");
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto result = conn.TableFunction("parquet_scan", {Value(path)})->Execute();
	timer.Note("parquet_scan opened");
	if (result->HasError()) {
		throw IOException("Failed to read parquet companion '%s': %s", path, result->GetError());
	}

	SampleInfo info;
	string label = "parquet companion '" + path + "'";
	IngestSampleResult(*result, label, info, /*load_iids=*/true, /*load_fids=*/true);
	info.sample_ct = info.iids.size();
	timer.Note("loaded %llu samples (iid_to_idx deferred)", (unsigned long long)info.sample_ct);
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
	BindPhaseTimer timer("LoadVariantMetadataFromSource");
	auto result = QuerySource(context, source, func_name);

	VariantMetadataIndex idx;
	idx.is_bim = false;
	string label = "source '" + source + "'";
	IngestVariantResult(*result, label, func_name, idx);
	idx.variant_ct = idx.chroms.size();
	idx.has_ids = true;
	idx.has_alleles = true;
	BuildChromOffsets(idx);
	return idx;
}

SampleInfo LoadSampleInfoFromSource(ClientContext &context, const string &source) {
	BindPhaseTimer timer("LoadSampleInfoFromSource");
	auto result = QuerySource(context, source, "LoadSampleInfoFromSource");
	SampleInfo info;
	string label = "source '" + source + "'";
	IngestSampleResult(*result, label, info, /*load_iids=*/true, /*load_fids=*/true);
	info.sample_ct = info.iids.size();
	return info;
}

// ---------------------------------------------------------------------------
// Unified dispatch functions
// ---------------------------------------------------------------------------

VariantMetadataIndex LoadVariantMetadata(ClientContext &context, const string &path, const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadata(dispatch:" + path + ")");
	if (IsParquetFile(path)) {
		return LoadVariantMetadataFromParquet(context, path, func_name);
	}
	if (IsNativePlinkFormat(path)) {
		return LoadVariantMetadataIndex(context, path, func_name);
	}
	return LoadVariantMetadataFromSource(context, path, func_name);
}

SampleInfo LoadSampleMetadata(ClientContext &context, const string &path) {
	BindPhaseTimer timer("LoadSampleMetadata(dispatch:" + path + ")");
	if (IsParquetFile(path)) {
		return LoadSampleInfoFromParquet(context, path);
	}
	if (IsNativePlinkFormat(path)) {
		return LoadSampleInfo(context, path);
	}
	return LoadSampleInfoFromSource(context, path);
}

SampleInfo LoadSampleCount(ClientContext &context, const string &path) {
	BindPhaseTimer timer("LoadSampleCount(" + path + ")");
	SampleInfo info;
	if (IsParquetFile(path)) {
		info.sample_ct = GetParquetRowCount(context, path);
		timer.Note("parquet metadata count = %llu", (unsigned long long)info.sample_ct);
		return info;
	}
	// Text/source: no efficient metadata-only count, fall through to a full load
	// then drop the IIDs to keep RSS predictable.
	auto full = IsNativePlinkFormat(path) ? LoadSampleInfo(context, path) : LoadSampleInfoFromSource(context, path);
	info.sample_ct = full.sample_ct;
	timer.Note("text fallback count = %llu (iids discarded)", (unsigned long long)info.sample_ct);
	return info;
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
		// Lazily build iid_to_idx — we only pay the ~500ms map-build cost at 7M
		// samples when VARCHAR sample filter is actually used.
		const_cast<SampleInfo *>(sample_info)->EnsureIidMap(func_name);
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

	VariantRange range;
	range.has_filter = true;

	// Fast path: chrom_offsets + binary search on POS (O(log N) per bound).
	auto it = variants.chrom_offsets.find(chrom);
	if (it == variants.chrom_offsets.end()) {
		return range; // empty
	}
	idx_t lo_local = it->second.first;
	idx_t hi_local = it->second.second;

	auto lb_local = [&](int32_t target) {
		idx_t lo = lo_local, hi = hi_local;
		while (lo < hi) {
			idx_t mid = lo + (hi - lo) / 2;
			if (variants.positions[mid] < target) {
				lo = mid + 1;
			} else {
				hi = mid;
			}
		}
		return lo;
	};
	idx_t start_local = lb_local(static_cast<int32_t>(start_pos));
	// Upper bound for end: first index with pos > end_pos.
	idx_t end_local = lo_local;
	{
		idx_t lo = lo_local, hi = hi_local;
		while (lo < hi) {
			idx_t mid = lo + (hi - lo) / 2;
			if (variants.positions[mid] <= static_cast<int32_t>(end_pos)) {
				lo = mid + 1;
			} else {
				hi = mid;
			}
		}
		end_local = lo;
	}

	// VariantRange is expressed in file-row vidx. In dense mode local == vidx.
	// Sparse indexes shouldn't hit this path (they've been pre-filtered).
	if (!variants.vidx_map.empty() && start_local < end_local) {
		throw InternalException("ParseRegion: sparse variant index should not re-filter by region");
	}
	range.start_idx = static_cast<uint32_t>(start_local);
	range.end_idx = static_cast<uint32_t>(end_local);
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

//! Lookup by (chrom, pos[, ref, alt]) using chrom_offsets + binary search on POS.
//! Returns the file-row vidx. Falls back to linear scan if chrom_offsets unavailable.
static uint32_t ResolveByCpra(const VariantMetadataIndex &variants, const string &chrom, int32_t pos,
                              const string *ref_match, const string *alt_match, const string &desc,
                              const string &func_name) {
	idx_t lo_local = 0;
	idx_t hi_local = variants.chroms.size();
	if (!variants.chrom_offsets.empty()) {
		auto it = variants.chrom_offsets.find(chrom);
		if (it == variants.chrom_offsets.end()) {
			throw InvalidInputException("%s: variant '%s' not found", func_name, desc);
		}
		lo_local = it->second.first;
		hi_local = it->second.second;
	}

	// Binary search on positions[lo_local..hi_local) — POS is sorted ascending within chrom.
	idx_t lo = lo_local, hi = hi_local;
	while (lo < hi) {
		idx_t mid = lo + (hi - lo) / 2;
		if (variants.positions[mid] < pos) {
			lo = mid + 1;
		} else {
			hi = mid;
		}
	}
	// Equal-range scan: there may be multiple variants at the same (chrom, pos) differing in alleles.
	for (idx_t i = lo; i < hi_local && variants.positions[i] == pos; i++) {
		if (variants.chrom_offsets.empty() && variants.chroms[i] != chrom) {
			continue; // fallback path: chrom check needed
		}
		if (ref_match && alt_match) {
			if (variants.refs[i] == *ref_match && variants.alts[i] == *alt_match) {
				// Map local idx back to file-row vidx
				if (!variants.vidx_map.empty()) {
					for (auto &kv : variants.vidx_map) {
						if (kv.second == i) {
							return kv.first;
						}
					}
					throw InternalException("ResolveByCpra: local idx %llu missing from vidx_map",
					                        static_cast<unsigned long long>(i));
				}
				return static_cast<uint32_t>(i);
			}
		} else {
			if (!variants.vidx_map.empty()) {
				for (auto &kv : variants.vidx_map) {
					if (kv.second == i) {
						return kv.first;
					}
				}
				throw InternalException("ResolveByCpra: local idx %llu missing from vidx_map",
				                        static_cast<unsigned long long>(i));
			}
			return static_cast<uint32_t>(i);
		}
	}

	throw InvalidInputException("%s: variant '%s' not found", func_name, desc);
}

//! Resolve a single CPRA string (chrom:pos or chrom:pos:ref:alt) to a variant index.
static uint32_t ResolveCpraString(const string &cpra, const VariantMetadataIndex &variants, const string &func_name) {
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

	char *end_ptr;
	errno = 0;
	long pos = std::strtol(parts[1].c_str(), &end_ptr, 10);
	if (*end_ptr != '\0' || errno != 0) {
		throw InvalidInputException("%s: invalid position in CPRA '%s'", func_name, cpra);
	}

	bool match_alleles = (parts.size() == 4);
	return ResolveByCpra(variants, parts[0], static_cast<int32_t>(pos), match_alleles ? &parts[2] : nullptr,
	                     match_alleles ? &parts[3] : nullptr, cpra, func_name);
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
	string desc = chrom + ":" + std::to_string(pos);
	if (match_alleles) {
		desc += ":" + ref_val + ":" + alt_val;
	}
	return ResolveByCpra(variants, chrom, pos, match_alleles ? &ref_val : nullptr, match_alleles ? &alt_val : nullptr,
	                     desc, func_name);
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
