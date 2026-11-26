#include "plink_common.hpp"

#include <atomic>
#include "duckdb_compat.hpp"
#include "plink_profile.hpp"

#include "duckdb/common/printer.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>

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

//! Build chrom_offsets from the chroms column. Requires (CHROM, POS)-sorted
//! input (PLINK spec). Detects non-contiguous chroms (a spec violation that
//! generic sources like SQL views might produce) and throws rather than
//! silently dropping runs. The rest of the code (binary search on POS) assumes
//! contiguous runs, so this is a hard precondition.
static void BuildChromOffsets(VariantMetadataIndex &idx, const string &source_label) {
	idx.chrom_offsets.clear();
	if (idx.chroms.empty()) {
		return;
	}
	idx_t run_start = 0;
	const string *current = &idx.chroms[0];
	for (idx_t i = 1; i < idx.chroms.size(); i++) {
		if (idx.chroms[i] != *current) {
			auto inserted = idx.chrom_offsets.emplace(*current, std::make_pair(run_start, i));
			if (!inserted.second) {
				throw InvalidInputException("%s: chromosome '%s' appears in non-contiguous runs (variants must be "
				                            "sorted by (CHROM, POS) as required by the PLINK spec; first run ended at "
				                            "row %llu, found again at row %llu)",
				                            source_label, current->c_str(),
				                            static_cast<unsigned long long>(inserted.first->second.second),
				                            static_cast<unsigned long long>(i));
			}
			current = &idx.chroms[i];
			run_start = i;
		}
	}
	auto inserted = idx.chrom_offsets.emplace(*current, std::make_pair(run_start, idx.chroms.size()));
	if (!inserted.second) {
		throw InvalidInputException(
		    "%s: chromosome '%s' appears in non-contiguous runs (variants must be sorted by (CHROM, POS))",
		    source_label, current->c_str());
	}
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
		// Track whether each required field was parsed (rather than left at default).
		// Prevents silent data corruption on truncated/malformed lines.
		bool got_chrom = false, got_pos = false, got_id = false, got_ref = false, got_alt = false;
		for (idx_t f = 0; f <= max_field; f++) {
			if (!NextField(buf, content_end, &cursor, &fstart, &flen, whitespace)) {
				break;
			}
			if (f == chrom_field) {
				chrom.assign(buf + fstart, flen);
				got_chrom = true;
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
				got_pos = true;
			} else if (f == id_field) {
				if (flen == 1 && buf[fstart] == '.') {
					// empty
				} else {
					id.assign(buf + fstart, flen);
				}
				got_id = true;
			} else if (f == ref_field) {
				ref.assign(buf + fstart, flen);
				got_ref = true;
			} else if (f == alt_field) {
				if (flen == 1 && buf[fstart] == '.') {
					// empty
				} else {
					alt.assign(buf + fstart, flen);
				}
				got_alt = true;
			}
		}

		// Reject truncated/malformed lines (missing required columns) instead of
		// silently pushing default-initialized values.
		if (!got_chrom || !got_pos || !got_id || !got_ref || !got_alt) {
			string missing;
			if (!got_chrom) {
				missing += " CHROM";
			}
			if (!got_pos) {
				missing += " POS";
			}
			if (!got_id) {
				missing += " ID";
			}
			if (!got_ref) {
				missing += " REF";
			}
			if (!got_alt) {
				missing += " ALT";
			}
			throw InvalidInputException(
			    "%s: .pvar/.bim file '%s' has a line missing required fields [%s] at byte offset %llu", func_name, path,
			    missing.c_str(), static_cast<unsigned long long>(pos));
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
	BuildChromOffsets(idx, path);
	return idx;
}

static void ParseVariantHeaderLine(const string &header_line, VariantMetadataIndex &idx, idx_t &chrom_field,
                                   idx_t &pos_field, idx_t &id_field, idx_t &ref_field, idx_t &alt_field) {
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
	} else {
		idx.is_bim = true;
		// .bim: CHROM(0) ID(1) CM(2) POS(3) ALT(4) REF(5)
		chrom_field = 0;
		id_field = 1;
		pos_field = 3;
		alt_field = 4;
		ref_field = 5;
	}
}

static bool ParseVariantDataLine(const char *buf, size_t content_end, bool whitespace, idx_t chrom_field,
                                 idx_t pos_field, idx_t id_field, idx_t ref_field, idx_t alt_field, const string &path,
                                 const string &func_name, uint64_t line_number, string &chrom, int32_t &pos_val,
                                 string &id, string &ref, string &alt) {
	const idx_t max_field = std::max({chrom_field, pos_field, id_field, ref_field, alt_field});
	size_t cursor = 0;
	size_t fstart = 0, flen = 0;
	bool got_chrom = false, got_pos = false, got_id = false, got_ref = false, got_alt = false;
	for (idx_t f = 0; f <= max_field; f++) {
		if (!NextField(buf, content_end, &cursor, &fstart, &flen, whitespace)) {
			break;
		}
		if (f == chrom_field) {
			chrom.assign(buf + fstart, flen);
			got_chrom = true;
		} else if (f == pos_field) {
			char tmp[32];
			idx_t n = flen < sizeof(tmp) - 1 ? flen : sizeof(tmp) - 1;
			std::memcpy(tmp, buf + fstart, n);
			tmp[n] = '\0';
			char *end;
			errno = 0;
			long v = std::strtol(tmp, &end, 10);
			if (end == tmp || *end != '\0' || errno != 0) {
				throw InvalidInputException("%s: invalid POS value '%s' at line %llu in '%s'", func_name, tmp,
				                            static_cast<unsigned long long>(line_number), path);
			}
			pos_val = static_cast<int32_t>(v);
			got_pos = true;
		} else if (f == id_field) {
			if (!(flen == 1 && buf[fstart] == '.')) {
				id.assign(buf + fstart, flen);
			}
			got_id = true;
		} else if (f == ref_field) {
			ref.assign(buf + fstart, flen);
			got_ref = true;
		} else if (f == alt_field) {
			if (!(flen == 1 && buf[fstart] == '.')) {
				alt.assign(buf + fstart, flen);
			}
			got_alt = true;
		}
	}

	if (!got_chrom || !got_pos || !got_id || !got_ref || !got_alt) {
		string missing;
		if (!got_chrom) {
			missing += " CHROM";
		}
		if (!got_pos) {
			missing += " POS";
		}
		if (!got_id) {
			missing += " ID";
		}
		if (!got_ref) {
			missing += " REF";
		}
		if (!got_alt) {
			missing += " ALT";
		}
		throw InvalidInputException("%s: .pvar/.bim file '%s' line %llu missing required fields [%s]", func_name, path,
		                            static_cast<unsigned long long>(line_number), missing.c_str());
	}
	return true;
}

VariantMetadataIndex LoadVariantMetadataFromTextRegion(ClientContext &context, const string &path, const string &chrom,
                                                       int64_t pos_start, int64_t pos_end, idx_t variant_ct_hint,
                                                       const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataFromTextRegion");
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	auto file_size = handle->GetFileSize();
	timer.Note("file_size=%llu", static_cast<unsigned long long>(file_size));

	if (file_size == 0) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is empty", func_name, path);
	}

	VariantMetadataIndex idx;
	// The vectors hold only the region-matching rows, indexed by local position
	// (file-row vidx via vidx_map / the contiguous fast path). Mark it a subset so
	// a zero-match region (empty vidx_map, no contiguous range) is not
	// misclassified as dense — see VariantMetadataIndex::is_subset. Mirrors
	// LoadVariantMetadataFromParquetRegion.
	idx.is_subset = true;
	idx.variant_ct = variant_ct_hint;
	idx.has_ids = true;
	idx.has_alleles = true;
	idx.chroms.reserve(1024);
	idx.positions.reserve(1024);
	idx.ids.reserve(1024);
	idx.refs.reserve(1024);
	idx.alts.reserve(1024);
	idx.local_to_vidx.reserve(1024);

	idx_t chrom_field = DConstants::INVALID_INDEX;
	idx_t pos_field = DConstants::INVALID_INDEX;
	idx_t id_field = DConstants::INVALID_INDEX;
	idx_t ref_field = DConstants::INVALID_INDEX;
	idx_t alt_field = DConstants::INVALID_INDEX;
	bool layout_ready = false;
	bool saw_target_chrom = false;
	uint32_t data_vidx = 0;
	uint64_t line_number = 0;

	// Stream fixed-size chunks and parse complete lines in-place. This keeps the
	// native text region loader bounded-memory while avoiding the extreme syscall
	// overhead of FileHandle::ReadLine() on 75M-row whole-genome PVARs.
	auto process_line = [&](const char *line_buf, size_t line_len) -> bool {
		line_number++;
		if (line_len > 0 && line_buf[line_len - 1] == '\r') {
			line_len--;
		}
		if (line_len == 0) {
			return true;
		}

		if (!layout_ready) {
			if (line_len >= 2 && line_buf[0] == '#' && line_buf[1] == '#') {
				return true;
			}
			string header_line(line_buf, line_len);
			ParseVariantHeaderLine(header_line, idx, chrom_field, pos_field, id_field, ref_field, alt_field);
			if (chrom_field == DConstants::INVALID_INDEX || pos_field == DConstants::INVALID_INDEX ||
			    id_field == DConstants::INVALID_INDEX || ref_field == DConstants::INVALID_INDEX ||
			    alt_field == DConstants::INVALID_INDEX) {
				throw InvalidInputException("%s: .pvar/.bim file '%s' is missing required columns "
				                            "(need CHROM, POS, ID, REF, ALT)",
				                            func_name, path);
			}
			layout_ready = true;
			if (!idx.is_bim) {
				return true; // header consumed; next line is data
			}
			// .bim has no header; fall through and parse this line as data.
		}

		if (!idx.is_bim && line_buf[0] == '#') {
			return true;
		}

		string row_chrom, row_id, row_ref, row_alt;
		int32_t row_pos = 0;
		ParseVariantDataLine(line_buf, line_len, idx.is_bim, chrom_field, pos_field, id_field, ref_field, alt_field,
		                     path, func_name, line_number, row_chrom, row_pos, row_id, row_ref, row_alt);

		bool chrom_match = (row_chrom == chrom);
		if (chrom_match) {
			saw_target_chrom = true;
			if (row_pos >= pos_start && row_pos <= pos_end) {
				uint32_t local = static_cast<uint32_t>(idx.chroms.size());
				idx.vidx_map.emplace(data_vidx, local);
				idx.local_to_vidx.push_back(data_vidx);
				idx.chroms.emplace_back(std::move(row_chrom));
				idx.positions.push_back(row_pos);
				idx.ids.emplace_back(std::move(row_id));
				idx.refs.emplace_back(std::move(row_ref));
				idx.alts.emplace_back(std::move(row_alt));
			} else if (row_pos > pos_end) {
				return false; // sorted PVAR/BIM: interval is exhausted
			}
		} else if (saw_target_chrom) {
			return false; // sorted by chromosome; target chromosome is exhausted
		}
		data_vidx++;
		return true;
	};

	static constexpr size_t READ_CHUNK_SIZE = 8ULL * 1024ULL * 1024ULL;
	vector<char> read_buffer(READ_CHUNK_SIZE);
	string carry;
	string line_buffer;
	line_buffer.reserve(READ_CHUNK_SIZE + 4096);
	uint64_t remaining = static_cast<uint64_t>(file_size);
	bool continue_scan = true;
	while (remaining > 0 && continue_scan) {
		size_t read_size = static_cast<size_t>(std::min<uint64_t>(remaining, READ_CHUNK_SIZE));
		handle->Read(read_buffer.data(), read_size);
		remaining -= read_size;

		line_buffer.clear();
		if (!carry.empty()) {
			line_buffer.append(carry);
		}
		line_buffer.append(read_buffer.data(), read_size);

		size_t line_start = 0;
		while (line_start < line_buffer.size()) {
			auto newline = line_buffer.find('\n', line_start);
			if (newline == string::npos) {
				break;
			}
			if (!process_line(line_buffer.data() + line_start, newline - line_start)) {
				continue_scan = false;
				break;
			}
			line_start = newline + 1;
		}
		if (continue_scan) {
			carry.assign(line_buffer.data() + line_start, line_buffer.size() - line_start);
		}
	}
	if (continue_scan && !carry.empty()) {
		process_line(carry.data(), carry.size());
	}

	if (!layout_ready) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' contains no header or data", func_name, path);
	}
	if (!idx.chroms.empty()) {
		idx.chrom_offsets.emplace(idx.chroms.front(), std::make_pair(idx_t {0}, idx.chroms.size()));
		if (!idx.local_to_vidx.empty() &&
		    idx.local_to_vidx.back() - idx.local_to_vidx.front() + 1 == idx.local_to_vidx.size()) {
			idx.has_contiguous_vidx_range = true;
			idx.contiguous_start_vidx = idx.local_to_vidx.front();
		}
	}
	timer.Note("ingested %llu region variants (variant_ct hint=%llu)", (unsigned long long)idx.chroms.size(),
	           (unsigned long long)variant_ct_hint);
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

vector<string> ResolvePathList(const Value &input, const char *fn_name) {
	vector<string> paths;
	// A NULL list value (e.g. list(x) over zero rows, or NULL::VARCHAR[]) has a LIST
	// type but no children — guard before GetChildren, which asserts on a NULL value.
	// (The child.IsNull() skip below only handles NULL *elements* of a non-null list.)
	if (input.IsNull()) {
		throw InvalidInputException("%s: empty file list provided", fn_name);
	}
	if (input.type().id() == LogicalTypeId::LIST) {
		for (auto &child : ListValue::GetChildren(input)) {
			if (!child.IsNull()) {
				paths.push_back(child.GetValue<string>());
			}
		}
		if (paths.empty()) {
			throw InvalidInputException("%s: empty file list provided", fn_name);
		}
	} else {
		paths.push_back(input.GetValue<string>());
	}
	return paths;
}

string ResolveExistingPath(ClientContext &context, FileSystem &fs, const string &path) {
	// Literal / cwd-relative / absolute path wins, matching DuckDB's precedence
	// in LocalFileSystem::FetchFileWithoutGlob.
	if (fs.FileExists(path)) {
		return path;
	}
	// Only relative paths get file_search_path expansion (matches DuckDB).
	if (fs.IsPathAbsolute(path)) {
		return "";
	}
	Value value;
	if (context.TryGetCurrentSetting("file_search_path", value)) {
		auto search_paths_str = value.ToString();
		if (!search_paths_str.empty()) {
			for (auto &search_path : StringUtil::Split(search_paths_str, ',')) {
				auto joined = fs.JoinPath(search_path, path);
				if (fs.FileExists(joined)) {
					return joined;
				}
			}
		}
	}
	return "";
}

vector<string> ExpandPathInputs(ClientContext &context, FileSystem &fs, const vector<string> &inputs,
                                const char *fn_name) {
	// NB: FileSystem::GetFileSystem(context) returns an OpenerFileSystem that
	// pushes the context-bound opener into Glob automatically — passing our own
	// opener is an error ("cannot take an opener"). So call Glob with no opener;
	// protocol filesystems that need the client (e.g. scalarfs `pathmacro:`,
	// which runs a catalog macro) still receive it via that pushed opener.
	(void)context;
	vector<string> expanded;
	for (auto &input : inputs) {
		if (input.empty()) {
			expanded.push_back(input); // named-params-only / empty positional
			continue;
		}
		const bool is_glob = FileSystem::HasGlob(input);
		// Route through the VFS so glob patterns and protocol URLs (e.g.
		// pathmacro:) resolve to concrete local paths, mirroring read_csv. Local
		// .pgen I/O then opens those concrete paths via pgenlib's raw FILE*.
		auto matches = fs.Glob(input);
		if (!matches.empty()) {
			// Glob order is filesystem-dependent (LocalFileSystem does not sort),
			// so sort each pattern's matches for a deterministic, reproducible
			// shard order. This matters for row order and especially by-index
			// variant selection, whose global offsets depend on shard order.
			// read_csv sorts its glob expansion for the same reason. Explicit
			// list order (across separate inputs) is preserved — only the
			// matches WITHIN one pattern are sorted.
			std::sort(matches.begin(), matches.end(),
			          [](const OpenFileInfo &a, const OpenFileInfo &b) { return a.path < b.path; });
			for (auto &m : matches) {
				expanded.push_back(m.path);
			}
		} else if (is_glob) {
			throw InvalidInputException("%s: no files matched pattern '%s'", fn_name, input);
		} else {
			// Plain prefix / literal path with no glob match: keep as-is so the
			// prefix->.pgen discovery and file_search_path resolution still run.
			expanded.push_back(input);
		}
	}
	return expanded;
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
//!
//! `expect_ordinal` selects how the per-variant ordinal is determined:
//!  * false (the dense full read / generic source) — the index is DENSE, numbered
//!    by emit order. A `file_row_number` column, if the source happens to carry
//!    one, is IGNORED: a stray user column must never be adopted as a .pgen seek
//!    ordinal.
//!  * true (the region pushdown) — the result MUST carry an ordinal column named
//!    `file_row_number` (the writer's `__vidx`, aliased); the index is SPARSE
//!    (idx.vidx_map populated). Every ordinal becomes a .pgen seek index and the
//!    column may come from an untrusted companion, so each value is validated
//!    against `vidx_bound` (the block's raw variant count, exclusive upper bound):
//!    reject < 0, >= vidx_bound, or a duplicate within the block. The bound is
//!    mandatory — there is no "disable validation" sentinel.
static void IngestVariantResult(QueryResult &result, const string &source_label, const string &func_name,
                                VariantMetadataIndex &idx, bool expect_ordinal, idx_t vidx_bound) {
	// Copied to strings rather than bound by reference: QueryResult::names is
	// vector<Identifier> on DuckDB v2.0, and the matching below is plain
	// lowercase string comparison.
	auto col_names = CompatResultNames(result);
	auto &col_types = CompatResultTypes(result);
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
		} else if (expect_ordinal && lower == "file_row_number") {
			// Only bind the ordinal when the caller expects one. In the dense path a column that
			// happens to be named file_row_number is left as data, never used as a seek index.
			rn_col = i;
		}
	}

	if (chrom_col == DConstants::INVALID_INDEX || pos_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: %s missing required columns (need CHROM and POS, found: %s)", func_name,
		                            source_label, StringUtil::Join(col_names, ", "));
	}

	// When the caller expects an ordinal (the region pushdown), it must be present — otherwise the
	// sparse vidx_map would be empty and every variant would mis-seek. Fail loud rather than silently
	// falling back to a dense/emit-order index for a path that needs real pgen indices.
	if (expect_ordinal && rn_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: %s has no ordinal column (expected 'file_row_number')", func_name,
		                            source_label);
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

	// The ordinal column is read below via reinterpret_cast<const int64_t*>, so its physical width
	// must be 64-bit. The region query CASTs __vidx to BIGINT, but a foreign companion read through a
	// different path could present it otherwise; guard here before the cast.
	if (rn_col != DConstants::INVALID_INDEX) {
		auto rn_type_id = col_types[rn_col].id();
		if (rn_type_id != LogicalTypeId::BIGINT && rn_type_id != LogicalTypeId::UBIGINT) {
			throw InvalidInputException("%s: %s ordinal column (__vidx/file_row_number) has unsupported type %s "
			                            "(expected BIGINT)",
			                            func_name, source_label, col_types[rn_col].ToString());
		}
	}

	unique_ptr<DataChunk> chunk;
	while ((chunk = result.Fetch()) != nullptr && chunk->size() > 0) {
		idx_t n = chunk->size();

		// CHROM — may be VARCHAR or numeric
		{
			auto &vec = chunk->data[chrom_col];
			UnifiedVectorFormat uvf;
			CompatToUnifiedFormat(vec, n, uvf);
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
			CompatToUnifiedFormat(vec, n, uvf);
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
				CompatToUnifiedFormat(vec, n, uvf);
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
				CompatToUnifiedFormat(vec, n, uvf);
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
				CompatToUnifiedFormat(vec, n, uvf);
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

		// file_row_number → sparse vidx_map + local_to_vidx reverse map.
		// Both are populated so callers iterating the loaded subset can map
		// local → vidx in O(1) (e.g. ResolveByCpra, BuildVariantIdIndex).
		if (rn_col != DConstants::INVALID_INDEX) {
			auto &vec = chunk->data[rn_col];
			UnifiedVectorFormat uvf;
			CompatToUnifiedFormat(vec, n, uvf);
			auto data = reinterpret_cast<const int64_t *>(uvf.data);
			idx_t local_base = idx.chroms.size() - n;
			if (idx.local_to_vidx.size() < idx.chroms.size()) {
				idx.local_to_vidx.resize(idx.chroms.size());
			}
			for (idx_t i = 0; i < n; i++) {
				auto si = uvf.sel->get_index(i);
				if (!uvf.validity.RowIsValid(si)) {
					throw IOException("%s: NULL file_row_number encountered (internal error)", func_name);
				}
				// This ordinal becomes a .pgen seek index and may come from an untrusted companion,
				// so bound it to [0, vidx_bound) before narrowing to uint32_t. An out-of-range value
				// would seek past the .pgen (OOB read); a duplicate would double-map a record. The bound
				// is enforced unconditionally — reached only when expect_ordinal is set, where vidx_bound
				// is the real raw variant count (an empty block has no rows here, so nothing to check).
				int64_t raw_vidx = data[si];
				if (raw_vidx < 0 || static_cast<uint64_t>(raw_vidx) >= vidx_bound) {
					throw InvalidInputException(
					    "%s: %s has __vidx %lld out of range [0, %llu) — companion does not match the .pgen", func_name,
					    source_label, static_cast<long long>(raw_vidx), static_cast<unsigned long long>(vidx_bound));
				}
				uint32_t vidx = static_cast<uint32_t>(raw_vidx);
				uint32_t local = static_cast<uint32_t>(local_base + i);
				auto inserted = idx.vidx_map.emplace(vidx, local);
				if (!inserted.second) {
					throw InvalidInputException("%s: %s has duplicate __vidx %u — companion is malformed", func_name,
					                            source_label, vidx);
				}
				idx.local_to_vidx[local] = vidx;
			}
		}
	}
}

VariantMetadataIndex LoadVariantMetadataFromParquet(ClientContext &context, const string &path,
                                                    const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataFromParquet");
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	// SQL, not Connection::TableFunction: DuckDB main removed the C++ relational API.
	auto escaped_path = StringUtil::Replace(path, "'", "''");
	auto result = conn.Query("SELECT * FROM parquet_scan('" + escaped_path + "')");
	timer.Note("parquet_scan opened");
	if (result->HasError()) {
		throw IOException("%s: failed to read parquet companion '%s': %s", func_name, path, result->GetError());
	}

	VariantMetadataIndex idx;
	idx.is_bim = false;
	string label = "parquet companion '" + path + "'";
	// Dense full read: numbered by emit order; any ordinal column is ignored, so no bound applies.
	IngestVariantResult(*result, label, func_name, idx, /*expect_ordinal=*/false, /*vidx_bound=*/0);
	idx.variant_ct = idx.chroms.size();
	idx.has_ids = true;
	idx.has_alleles = true;
	timer.Note("ingested %llu variants", (unsigned long long)idx.variant_ct);
	BuildChromOffsets(idx, path);
	timer.Note("built chrom_offsets (%llu chroms)", (unsigned long long)idx.chrom_offsets.size());
	return idx;
}

//! Pick the chromosome column name. Different pvar producers use either
//! "#CHROM" (plink2 native) or "CHROM" (custom exports). Introspecting via
//! a LIMIT 0 scan is cheaper than a wrapped DESCRIBE: DuckDB only reads the
//! parquet footer to get schema, no row groups touched.
//!
//! `has_vidx_out` reports whether the companion carries an explicit `__vidx`
//! column (the writer's 0-based pgen record index). When present it is the
//! authoritative ordinal source; the caller reads it instead of the physical
//! file_row_number, which is only the pgen index when the file is stored in
//! __vidx order (which a re-sorted or re-written companion need not be).
static string DetectChromColumn(Connection &conn, const string &path, const string &func_name, bool &has_vidx_out) {
	auto escaped_path = StringUtil::Replace(path, "'", "''");
	auto result = conn.Query("SELECT * FROM read_parquet('" + escaped_path + "') LIMIT 0");
	if (result->HasError()) {
		throw IOException("%s: failed to read parquet schema for '%s': %s", func_name, path, result->GetError());
	}
	// BaseQueryResult::names is private on v2.0; reach it through the shim.
	auto schema_names = CompatResultNames(*result);
	has_vidx_out = false;
	for (auto &name : schema_names) {
		if (StringUtil::Lower(name) == "__vidx") {
			has_vidx_out = true;
			break;
		}
	}
	for (auto &name : schema_names) {
		auto lower = StringUtil::Lower(name);
		if (lower == "#chrom") {
			return "#CHROM";
		}
	}
	for (auto &name : schema_names) {
		auto lower = StringUtil::Lower(name);
		if (lower == "chrom") {
			return "CHROM";
		}
	}
	throw InvalidInputException("%s: parquet companion '%s' has no CHROM or #CHROM column (columns: %s)", func_name,
	                            path, StringUtil::Join(schema_names, ", "));
}

//! Region-pushdown loader: queries the parquet file with a WHERE clause so
//! only matching rows are materialized. Returns a sparse index whose
//! vidx_map keys are the file row numbers (pgenlib-compatible vidx).
//!
//! `variant_ct_hint` — the caller's expected total row count (typically
//! `pgfi.raw_variant_ct` from the pgen header). Used as `idx.variant_ct` so
//! the downstream count-mismatch check is effectively a sanity check against
//! pgen (we don't re-scan the full parquet to count rows; trusting pgen keeps
//! the region path O(region size), not O(file size)).
VariantMetadataIndex LoadVariantMetadataFromParquetRegion(ClientContext &context, const string &path,
                                                          const string &chrom, int64_t pos_start, int64_t pos_end,
                                                          idx_t variant_ct_hint, const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataFromParquetRegion");
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Each query on the companion opens it again, and on a remote file set every open is a HEAD and
	// a round trip (#127). So the loader runs ONE query, for the layout the plinking writer always
	// produces: a "CHROM" column and an explicit __vidx. Only when that query fails does it read the
	// schema (DetectChromColumn) and choose the layout. If the schema IS the assumed layout, the first
	// error is a real one (POS type, a bad file) and is thrown unchanged, so no schema error is masked
	// as a column-name failure. One plain query, not Prepare + Execute: each binds, and each bind
	// opens the file.
	//
	// One known difference from the probe-first loader: a file that has BOTH a "#CHROM" and a "CHROM"
	// column AND __vidx is filtered on "CHROM" here, where DetectChromColumn picks "#CHROM". The
	// "CHROM" column can also come from the path: a file under a CHROM=<v>/ directory gets one from
	// hive auto-detection. The plinking writer writes "CHROM" only, so this needs a converted file.
	// Without __vidx the first query fails and DetectChromColumn's choice is used, as before.
	//
	// The ordinal (pgenlib-compatible global vidx) must be the variant's index in
	// the .pgen, independent of the companion's physical row order, because row-group
	// pruning and a re-sorted companion both move rows. When the companion carries the
	// writer's explicit __vidx column, read THAT (aliased to file_row_number so the
	// ingest helper's by-name match is unchanged). Only fall back to the physical
	// file_row_number for a legacy companion that predates __vidx; there the pgen index
	// is recoverable only from physical position, which is correct for an in-order file.
	// Either way stats pruning + column projection keep this O(region_size), not O(total).
	//
	// The path and the user-controlled chrom are embedded via quote-doubling (DuckDB follows
	// SQL-standard string literals — no backslash escapes — so `''` is sufficient). The positions
	// are int64 values, printed as integers.
	auto escaped_path = StringUtil::Replace(path, "'", "''");
	auto escaped_chrom = StringUtil::Replace(chrom, "'", "''");
	auto region_sql = [&](const string &chrom_col_name, bool has_vidx) {
		// CAST __vidx to BIGINT: a no-op when the writer's companion already stores it as BIGINT, but it
		// normalizes a foreign companion that wrote __vidx as INT32/UINT so the int64 ingest read is sound.
		string ordinal_expr = has_vidx ? "CAST(__vidx AS BIGINT) AS file_row_number" : "file_row_number";
		string read_parquet_args = has_vidx ? "" : ", file_row_number=true";
		string order_by = has_vidx ? "__vidx" : "file_row_number";
		return "SELECT \"" + chrom_col_name + "\" AS \"#CHROM\", POS, ID, REF, ALT, " + ordinal_expr +
		       " FROM read_parquet('" + escaped_path + "'" + read_parquet_args + ") WHERE (CAST(\"" + chrom_col_name +
		       "\" AS VARCHAR) = '" + escaped_chrom + "') AND POS BETWEEN " + std::to_string(pos_start) + " AND " +
		       std::to_string(pos_end) + " ORDER BY " + order_by;
	};
	const string writer_chrom_col = "CHROM";
	auto result = conn.Query(region_sql(writer_chrom_col, /*has_vidx=*/true));
	if (result->HasError()) {
		bool has_vidx = false;
		auto chrom_col_name = DetectChromColumn(conn, path, func_name, has_vidx);
		timer.Note("chrom column: %s; __vidx: %s", chrom_col_name.c_str(), has_vidx ? "yes" : "no");
		if (chrom_col_name == writer_chrom_col && has_vidx) {
			throw IOException("%s: region-pushdown query failed on '%s': %s", func_name, path, result->GetError());
		}
		result = conn.Query(region_sql(chrom_col_name, has_vidx));
		if (result->HasError()) {
			throw IOException("%s: region-pushdown query failed on '%s': %s", func_name, path, result->GetError());
		}
	}
	timer.Note("pushdown query returned");

	VariantMetadataIndex idx;
	idx.is_bim = false;
	// The vectors hold only the region-matching rows, indexed by local position
	// (file-row vidx via vidx_map). Mark it a subset so a zero-match region (empty
	// vidx_map) is not misclassified as dense — see VariantMetadataIndex::is_subset.
	idx.is_subset = true;
	string label = "parquet companion '" + path + "' (region pushdown)";
	// Every ingested ordinal (__vidx, or a legacy file_row_number) must index within the .pgen's
	// variant_ct_hint variants; IngestVariantResult rejects an out-of-range or duplicate value.
	IngestVariantResult(*result, label, func_name, idx, /*expect_ordinal=*/true, /*vidx_bound=*/variant_ct_hint);
	// variant_ct = caller's hint (typically pgen's raw_variant_ct). The post-load
	// count-mismatch check in the caller becomes a sanity check against pgen;
	// we don't open the parquet file a second time to count rows.
	idx.variant_ct = variant_ct_hint;
	idx.has_ids = true;
	idx.has_alleles = true;
	// Single chrom in pushdown, single contiguous range of local indices
	if (!idx.chroms.empty()) {
		idx.chrom_offsets.emplace(idx.chroms.front(), std::make_pair(idx_t {0}, idx.chroms.size()));
		if (!idx.local_to_vidx.empty() &&
		    idx.local_to_vidx.back() - idx.local_to_vidx.front() + 1 == idx.local_to_vidx.size()) {
			idx.has_contiguous_vidx_range = true;
			idx.contiguous_start_vidx = idx.local_to_vidx.front();
		}
	}
	timer.Note("ingested %llu region variants (variant_ct hint=%llu)", (unsigned long long)idx.chroms.size(),
	           (unsigned long long)variant_ct_hint);
	return idx;
}

VariantMetadataIndex LoadVariantMetadataFromParquetOrdinals(ClientContext &context, const string &path,
                                                            const vector<uint32_t> &ordinals, idx_t variant_ct_hint,
                                                            const string &func_name) {
	BindPhaseTimer timer("LoadVariantMetadataFromParquetOrdinals");
	VariantMetadataIndex idx;
	idx.is_bim = false;
	idx.is_subset = true; // see LoadVariantMetadataFromParquetRegion
	idx.variant_ct = variant_ct_hint;
	idx.has_ids = true;
	idx.has_alleles = true;
	if (ordinals.empty()) {
		return idx;
	}
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	bool has_vidx = false;
	auto chrom_col_name = DetectChromColumn(conn, path, func_name, has_vidx);
	// The ordinals are integers, so the IN list is built as text; the path is quote-doubled.
	string in_list;
	for (auto o : ordinals) {
		in_list += (in_list.empty() ? "" : ", ") + std::to_string(o);
	}
	auto escaped_path = StringUtil::Replace(path, "'", "''");
	const string ordinal_col = has_vidx ? "CAST(__vidx AS BIGINT)" : "file_row_number";
	string sql = "SELECT \"" + chrom_col_name + "\" AS \"#CHROM\", POS, ID, REF, ALT, " + ordinal_col +
	             " AS file_row_number FROM read_parquet('" + escaped_path + "'" +
	             (has_vidx ? "" : ", file_row_number=true") + ") WHERE " + (has_vidx ? "__vidx" : "file_row_number") +
	             " IN (" + in_list + ") ORDER BY 6";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw IOException("%s: ordinal query failed on '%s': %s", func_name, path, result->GetError());
	}
	string label = "parquet companion '" + path + "' (ID index)";
	IngestVariantResult(*result, label, func_name, idx, /*expect_ordinal=*/true, /*vidx_bound=*/variant_ct_hint);
	BuildChromOffsets(idx, label);
	if (!idx.local_to_vidx.empty() &&
	    idx.local_to_vidx.back() - idx.local_to_vidx.front() + 1 == idx.local_to_vidx.size()) {
		idx.has_contiguous_vidx_range = true;
		idx.contiguous_start_vidx = idx.local_to_vidx.front();
	}
	timer.Note("ingested %llu of %llu requested ordinals", (unsigned long long)idx.chroms.size(),
	           (unsigned long long)ordinals.size());
	return idx;
}

//! Get total row count from a parquet file.
//!
//! Uses `SELECT COUNT(*) FROM read_parquet(...)`, which DuckDB optimizes to
//! read row-group `num_rows` from the footer — O(num_row_groups), not a
//! full scan. At biobank scale this is sub-millisecond (7M rows completed
//! in <1ms in our profiling), but it is NOT literally O(1): files with
//! many tiny row groups pay proportionally more. Use only for paths where
//! the whole file hasn't already been opened (e.g. count-only psam mode).
//! The region-pushdown path explicitly avoids this by trusting pgen's
//! raw_variant_ct as the variant_ct.
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

uint8_t CodeSampleSex(const Value &v) {
	if (v.IsNull()) {
		return 0;
	}
	try {
		int64_t s = v.DefaultCastAs(LogicalType::BIGINT).GetValue<int64_t>();
		return s == 1 ? 1 : (s == 2 ? 2 : 0);
	} catch (...) {
		return 0;
	}
}

//! SEX for ploidy-aware chrX/Y/MT stats. The column may be INTEGER (PLINK convention) or VARCHAR.
//! 1=male, 2=female, else 0: CodeSampleSex's rule, one vector cast per chunk (a Value per row took
//! most of 6.4 s for 8M samples, measured 2026-11-23). With an error string the try-cast sets a
//! failed row NULL instead of throwing, as CodeSampleSex's catch codes it 0.
static void AppendCodedSexes(Vector &sex_col, idx_t n, vector<uint8_t> &out) {
	Vector big(LogicalType::BIGINT, n);
	string cast_error;
	VectorOperations::DefaultTryCast(sex_col, big, n, &cast_error);
	UnifiedVectorFormat uvf;
	CompatToUnifiedFormat(big, n, uvf);
	auto data = UnifiedVectorFormat::GetData<int64_t>(uvf);
	for (idx_t i = 0; i < n; i++) {
		auto si = uvf.sel->get_index(i);
		const int64_t s = uvf.validity.RowIsValid(si) ? data[si] : 0;
		out.push_back(s == 1 ? 1 : (s == 2 ? 2 : 0));
	}
}

//! Columnar ingest helper for psam query results.
static void IngestSampleResult(QueryResult &result, const string &source_label, SampleInfo &info, bool load_iids,
                               bool load_fids) {
	// See IngestVariantResult: QueryResult::names is vector<Identifier> on v2.0.
	auto col_names = CompatResultNames(result);
	idx_t iid_col = DConstants::INVALID_INDEX;
	idx_t fid_col = DConstants::INVALID_INDEX;
	idx_t sex_col = DConstants::INVALID_INDEX;
	for (idx_t i = 0; i < col_names.size(); i++) {
		auto lower = StringUtil::Lower(col_names[i]);
		if (lower == "iid") {
			iid_col = i;
		} else if (lower == "fid") {
			fid_col = i;
		} else if (lower == "sex") {
			sex_col = i;
		}
	}
	if (iid_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s missing required IID column (found: %s)", source_label,
		                            StringUtil::Join(col_names, ", "));
	}
	bool has_fid_col = (fid_col != DConstants::INVALID_INDEX);
	bool has_sex_col = (sex_col != DConstants::INVALID_INDEX);

	idx_t total_rows = 0;
	unique_ptr<DataChunk> chunk;
	while ((chunk = result.Fetch()) != nullptr && chunk->size() > 0) {
		idx_t n = chunk->size();
		total_rows += n;

		// SEX for ploidy-aware chrX/Y/MT stats (see AppendCodedSexes).
		if (has_sex_col) {
			AppendCodedSexes(chunk->data[sex_col], n, info.sexes);
		}

		if (load_iids) {
			// An integer IID (a foreign companion: pandas writes int64) is read as its text, as a
			// text .psam reads it. Reading it as string_t crashed (SIGSEGV, 2026-11-23).
			Vector varchar_vec(LogicalType::VARCHAR, n);
			Vector *src = &chunk->data[iid_col];
			if (src->GetType().id() != LogicalTypeId::VARCHAR) {
				VectorOperations::DefaultCast(*src, varchar_vec, n);
				src = &varchar_vec;
			}
			auto &vec = *src;
			UnifiedVectorFormat uvf;
			CompatToUnifiedFormat(vec, n, uvf);
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
			// An integer FID (a foreign companion: pandas writes int64) is read as its text, as a
			// text .psam reads it. Reading it as string_t crashed (SIGSEGV, 2026-11-23).
			Vector varchar_vec(LogicalType::VARCHAR, n);
			Vector *src = &chunk->data[fid_col];
			if (src->GetType().id() != LogicalTypeId::VARCHAR) {
				VectorOperations::DefaultCast(*src, varchar_vec, n);
				src = &varchar_vec;
			}
			auto &vec = *src;
			UnifiedVectorFormat uvf;
			CompatToUnifiedFormat(vec, n, uvf);
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

static std::atomic<int64_t> g_psam_parquet_iid_loads {0};

void NotePsamParquetIidLoad() {
	g_psam_parquet_iid_loads.fetch_add(1, std::memory_order_relaxed);
}

int64_t PsamParquetIidLoads(bool reset) {
	if (reset) {
		g_psam_parquet_iid_loads.store(0, std::memory_order_relaxed);
		return 0;
	}
	return g_psam_parquet_iid_loads.load(std::memory_order_relaxed);
}

SampleInfo LoadSampleInfoFromParquet(ClientContext &context, const string &path) {
	BindPhaseTimer timer("LoadSampleInfoFromParquet");
	NotePsamParquetIidLoad();
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	// Read only the columns IngestSampleResult uses (IID, FID, SEX): the materialized
	// SELECT * held every column of every row (8M samples, measured 2026-11-23). Without an IID
	// column the scan keeps every column, so the "missing required IID column" error lists them.
	const string from = " FROM read_parquet('" + StringUtil::Replace(path, "'", "''") + "')";
	string select_list = "*";
	{
		auto probe = conn.Query("SELECT *" + from + " LIMIT 0");
		if (probe->HasError()) {
			throw IOException("Failed to read parquet companion '%s': %s", path, probe->GetError());
		}
		vector<string> keep;
		bool has_iid = false;
		for (auto &name : CompatResultNames(*probe)) {
			auto lower = StringUtil::Lower(name);
			if (lower == "iid" || lower == "fid" || lower == "sex") {
				has_iid |= lower == "iid";
				keep.push_back("\"" + StringUtil::Replace(name, "\"", "\"\"") + "\"");
			}
		}
		if (has_iid) {
			select_list = StringUtil::Join(keep, ", ");
		}
	}
	auto result = conn.Query("SELECT " + select_list + from);
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

bool LookupParquetIidRows(ClientContext &context, const string &path, const vector<string> &names,
                          vector<std::pair<string, uint32_t>> &matches, vector<uint8_t> *sexes, bool *has_sex) {
	if (names.empty() || names.size() > kMaxIidLookupNames || !IsParquetFile(path) || FileSystem::HasGlob(path)) {
		return false;
	}
	auto lit = [](const string &v) {
		return "'" + StringUtil::Replace(v, "'", "''") + "'";
	};
	auto quote = [](const string &v) {
		return "\"" + StringUtil::Replace(v, "\"", "\"\"") + "\"";
	};
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const string from = " FROM read_parquet(" + lit(path) + ", file_row_number := true)";
	auto probe = conn.Query("SELECT *" + from + " LIMIT 0");
	if (probe->HasError()) {
		return false; // the full load raises its own error
	}
	const auto col_names = CompatResultNames(*probe);
	const auto &col_types = CompatResultTypes(*probe);
	string iid_col;
	string sex_col;
	for (idx_t i = 0; i < col_names.size(); i++) {
		const auto lower = StringUtil::Lower(col_names[i]);
		if (lower == "iid") {
			if (col_types[i].id() != LogicalTypeId::VARCHAR) {
				return false;
			}
			iid_col = quote(col_names[i]);
		} else if (lower == "sex") {
			sex_col = quote(col_names[i]);
		}
	}
	if (iid_col.empty()) {
		return false;
	}
	const bool want_sex = sexes && !sex_col.empty();
	string in_list;
	for (idx_t k = 0; k < names.size(); k++) {
		in_list += (k ? ", " : "") + lit(names[k]);
	}
	auto result = conn.Query("SELECT file_row_number, " + iid_col + (want_sex ? ", " + sex_col : string()) + from +
	                         " WHERE " + iid_col + " IN (" + in_list + ") ORDER BY file_row_number");
	if (result->HasError()) {
		return false;
	}
	matches.clear();
	if (sexes) {
		sexes->clear();
	}
	if (has_sex) {
		*has_sex = want_sex;
	}
	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		for (idx_t r = 0; r < chunk->size(); r++) {
			matches.emplace_back(chunk->GetValue(1, r).GetValue<string>(),
			                     static_cast<uint32_t>(chunk->GetValue(0, r).GetValue<int64_t>()));
		}
		if (want_sex) {
			AppendCodedSexes(chunk->data[2], chunk->size(), *sexes);
		}
	}
	return true;
}

bool LookupNamedSamplesInParquet(ClientContext &context, const string &path, const Value &samples, SampleInfo &info,
                                 vector<std::pair<string, uint32_t>> &matches) {
	if (samples.IsNull() || samples.type().id() != LogicalTypeId::LIST ||
	    ListType::GetChildType(samples.type()).id() != LogicalTypeId::VARCHAR) {
		return false;
	}
	vector<string> names;
	for (auto &child : ListValue::GetChildren(samples)) {
		if (child.IsNull()) {
			return false;
		}
		names.push_back(child.GetValue<string>());
	}
	vector<uint8_t> sexes;
	bool has_sex = false;
	if (!LookupParquetIidRows(context, path, names, matches, &sexes, &has_sex)) {
		return false;
	}
	info = SampleInfo();
	info.sample_ct = GetParquetRowCount(context, path);
	if (has_sex) {
		info.sexes.assign(info.sample_ct, 0);
		for (idx_t m = 0; m < matches.size(); m++) {
			if (matches[m].second < info.sample_ct) {
				info.sexes[matches[m].second] = sexes[m];
			}
		}
	}
	return true;
}

Value NamedSampleIndices(const Value &samples, const vector<std::pair<string, uint32_t>> &matches,
                         const string &func_name) {
	unordered_map<string, vector<uint32_t>> rows;
	for (auto &m : matches) {
		rows[m.first].push_back(m.second);
	}
	vector<Value> indices;
	for (auto &child : ListValue::GetChildren(samples)) {
		auto name = child.GetValue<string>();
		auto it = rows.find(name);
		if (it == rows.end()) {
			throw InvalidInputException("%s: sample '%s' not found in .psam", func_name, name);
		}
		if (it->second.size() > 1) {
			throw IOException("%s has duplicate IID '%s' (at sample %llu and %llu)", func_name, name,
			                  static_cast<unsigned long long>(it->second[0] + 1),
			                  static_cast<unsigned long long>(it->second[1] + 1));
		}
		indices.push_back(Value::BIGINT(it->second[0]));
	}
	return Value::LIST(LogicalType::BIGINT, std::move(indices));
}

// ---------------------------------------------------------------------------
// Generic source loading (CSV, tables, views, etc.)
// ---------------------------------------------------------------------------

//! Execute a query against a source string via a separate Connection.
//! The source may be a file path (CSV, parquet, etc.) or a table/view name.
//! DuckDB's replacement scan mechanism handles file extension dispatch.
static unique_ptr<QueryResult> QuerySource(ClientContext &context, const string &source, const string &func_name) {
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
	// Generic source: dense, emit-order; a stray file_row_number column is not adopted as an ordinal.
	IngestVariantResult(*result, label, func_name, idx, /*expect_ordinal=*/false, /*vidx_bound=*/0);
	idx.variant_ct = idx.chroms.size();
	idx.has_ids = true;
	idx.has_alleles = true;
	BuildChromOffsets(idx, source);
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
	// NULL value (bare NULL, or a LIST-typed NULL from list(x) over zero rows) has no
	// children/child-type to read — guard before ListValue/ListType access.
	if (samples_val.IsNull()) {
		throw InvalidInputException("%s: samples list must not be empty", func_name);
	}
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
		if (!sample_info->iid_to_idx.empty()) {
			// A caller already built the map (read_pfile): use it.
			for (auto &child : children) {
				auto iid = child.GetValue<string>();
				auto it = sample_info->iid_to_idx.find(iid);
				if (it == sample_info->iid_to_idx.end()) {
					throw InvalidInputException("%s: sample '%s' not found in .psam", func_name, iid);
				}
				indices.push_back(static_cast<uint32_t>(it->second));
			}
		} else {
			if (sample_info->iids.empty() && sample_info->sample_ct > 0) {
				throw InternalException("%s: IID lookup requested but IID strings were not loaded (%llu samples)",
				                        func_name, static_cast<unsigned long long>(sample_info->sample_ct));
			}
			// One pass over the IIDs for the named ones, not a map of every IID: 8M samples,
			// plink_freq(samples := ['per5']) took 5.80 s with the map, 0.60 s by index (measured
			// 2026-11-24). Only a named IID is checked for a duplicate (as read_pfile since #150).
			unordered_map<string, vector<idx_t>> rows;
			for (auto &child : children) {
				rows.emplace(child.GetValue<string>(), vector<idx_t>());
			}
			const auto &iids = sample_info->iids;
			for (idx_t i = 0; i < iids.size(); i++) {
				auto it = rows.find(iids[i]);
				if (it != rows.end() && it->second.size() < 2) {
					it->second.push_back(i);
				}
			}
			for (auto &child : children) {
				auto iid = child.GetValue<string>();
				auto &found = rows[iid];
				if (found.empty()) {
					throw InvalidInputException("%s: sample '%s' not found in .psam", func_name, iid);
				}
				if (found.size() > 1) {
					throw IOException("%s has duplicate IID '%s' (at sample %llu and %llu)", func_name, iid,
					                  static_cast<unsigned long long>(found[0] + 1),
					                  static_cast<unsigned long long>(found[1] + 1));
				}
				indices.push_back(static_cast<uint32_t>(found[0]));
			}
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
                             const string &func_name, bool *include_missing_out) {
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

		if (include_missing_out && field_name == "include_missing") {
			if (!child_val.IsNull()) {
				*include_missing_out = child_val.GetValue<bool>();
			}
			continue;
		}

		if (field_name != "min" && field_name != "max") {
			throw InvalidInputException(
			    "%s: %s has unknown field '%s' (expected %s)", func_name, param_name, CompatNameStr(field_name),
			    include_missing_out ? "'min', 'max', and/or 'include_missing'" : "'min' and/or 'max'");
		}

		if (child_val.IsNull()) {
			continue; // NULL = unbounded
		}

		double v = child_val.GetValue<double>();
		if (v < valid_min || v > valid_max) {
			throw InvalidInputException("%s: %s.%s value %g is out of range [%g, %g]", func_name, param_name,
			                            CompatNameStr(field_name), v, valid_min, valid_max);
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

void ParseIncludeGenotypes(const Value &val, GenotypeRangeFilter &out, const string &func_name) {
	// NULL (bare, or a LIST-typed NULL) = no filter, matching the empty-list case below.
	if (val.IsNull()) {
		return;
	}
	if (val.type().id() != LogicalTypeId::LIST) {
		throw InvalidInputException("%s: include_genotypes must be a LIST of category names "
		                            "(e.g. ['het', 'hom_alt'])",
		                            func_name);
	}

	auto &children = ListValue::GetChildren(val);
	if (children.empty()) {
		return; // empty list = no filter
	}

	for (auto &child : children) {
		if (child.IsNull()) {
			throw InvalidInputException("%s: include_genotypes contains a NULL category name", func_name);
		}
		string label = child.GetValue<string>();
		StringUtil::Trim(label);
		label = StringUtil::Lower(label);
		if (label == "hom_ref") {
			out.allowed[0] = true;
		} else if (label == "het") {
			out.allowed[1] = true;
		} else if (label == "hom_alt") {
			out.allowed[2] = true;
		} else if (label == "missing") {
			out.include_missing = true;
		} else {
			throw InvalidInputException("%s: include_genotypes has unknown category '%s' "
			                            "(expected 'hom_ref', 'het', 'hom_alt', and/or 'missing')",
			                            func_name, label);
		}
	}

	out.active = true;
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

GenotypeRangeResult CheckGenotypeRange(const GenotypeRangeFilter &filter, const STD_ARRAY_REF(uint32_t, 4) genocounts) {
	GenotypeRangeResult result;
	result.any_pass = false;
	result.all_pass = true;

	// genocounts: [hom_ref(0), het(1), hom_alt(2), missing(3)]
	// Genotype values map directly to array indices 0, 1, 2
	for (int g = 0; g <= 2; g++) {
		bool selected = filter.allowed[g];
		if (selected && genocounts[g] > 0) {
			result.any_pass = true;
		}
		if (!selected && genocounts[g] > 0) {
			result.all_pass = false;
		}
	}

	// Missing samples are retained only when include_missing is set; when they are,
	// a variant whose only kept samples are missing must not be skipped by any_pass.
	// all_pass is unaffected: missing hardcalls are already emitted as NULL, so the
	// per-element null-out fast path does not need to account for them.
	if (filter.include_missing && genocounts[3] > 0) {
		result.any_pass = true;
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
		auto gr = CheckGenotypeRange(genotype_filter, genocounts);
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

unordered_map<string, uint32_t> BuildVariantIdIndex(const VariantMetadataIndex &variants) {
	// Iterate the loaded subset (dense or sparse) by local index. In sparse mode
	// the local slot's file-row vidx comes from local_to_vidx; in dense mode
	// local == vidx. Using variant_ct to bound the loop is WRONG for sparse
	// indexes (where variant_ct is the full source count, not the loaded size).
	unordered_map<string, uint32_t> id_to_idx;
	id_to_idx.reserve(variants.chroms.size());
	for (idx_t local = 0; local < variants.chroms.size(); local++) {
		const auto &id = variants.ids[local];
		if (!id.empty()) {
			id_to_idx[id] = variants.VidxForLocal(local);
		}
	}
	return id_to_idx;
}

uint64_t CountMultiallelicInRange(const VariantMetadataIndex &variants, uint32_t start_vidx, uint32_t end_vidx) {
	// A lazy parquet path may have skipped REF/ALT entirely; alts is then empty
	// (or stale) and there is nothing to inspect.
	if (!variants.has_alleles) {
		return 0;
	}
	// Iterate the LOADED subset, exactly as BuildVariantIdIndex does. Walking
	// vidx from start_vidx to end_vidx instead would call Local() on rows that a
	// region pushdown never loaded, which throws. chroms is the canonical subset
	// length; alts is only indexed where it actually has a row.
	uint64_t multiallelic_ct = 0;
	const auto subset_size = MinValue<idx_t>(variants.chroms.size(), variants.alts.size());
	for (idx_t local = 0; local < subset_size; local++) {
		const auto vidx = variants.VidxForLocal(local);
		if (vidx < start_vidx || vidx >= end_vidx) {
			continue;
		}
		if (variants.alts[local].find(',') != string::npos) {
			multiallelic_ct++;
		}
	}
	return multiallelic_ct;
}

void WarnOnMultiallelicCollapse(const VariantMetadataIndex &variants, uint32_t start_vidx, uint32_t end_vidx,
                                const char *func_name) {
	WarnOnMultiallelicCount(CountMultiallelicInRange(variants, start_vidx, end_vidx), func_name);
}

void WarnOnMultiallelicCount(uint64_t multiallelic_ct, const char *func_name) {
	if (multiallelic_ct == 0) {
		return;
	}
	Printer::Print(StringUtil::Format("%s: %llu variant(s) have multiple ALT alleles; statistics collapse them into "
	                                  "a single REF-vs-any-ALT figure and ALT is reported verbatim",
	                                  func_name, static_cast<unsigned long long>(multiallelic_ct)));
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
				return variants.VidxForLocal(i);
			}
		} else {
			return variants.VidxForLocal(i);
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
	// NULL value (bare NULL, or a LIST-typed NULL from list(x) over zero rows) has no
	// children/child-type — guard before the type dispatch touches ListValue/ListType.
	if (val.IsNull()) {
		throw InvalidInputException("%s: variants must not be NULL", func_name);
	}
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
// Unified axis selection handler (variant × sample axes)
// ---------------------------------------------------------------------------

//! Destructure a `selection := <STRUCT>` value's fields into a Selection. PURE:
//! reads the struct, no IO, no ordinal resolution. Fields are matched by EXACT name
//! (matching the struct's own field names, as ResolveCpraStruct does) so a typo
//! like `sample` is rejected as unknown rather than silently dropped. A NULL field
//! is treated as "unset" for that axis. These seed the lowest precedence; the
//! individual named params override them afterwards.
static void MergeSelectionStruct(const Value &select_val, Selection &sel) {
	if (select_val.type().id() != LogicalTypeId::STRUCT) {
		throw InvalidInputException("selection: must be a STRUCT of selection fields "
		                            "(samples, variants, region, variant_range), got %s",
		                            select_val.type().ToString());
	}
	auto &child_types = StructType::GetChildTypes(select_val.type());
	auto &children = StructValue::GetChildren(select_val);
	for (idx_t i = 0; i < child_types.size(); i++) {
		const auto &field = child_types[i].first; // Identifier on v2.0, string on v1.5; compared below
		const Value &fv = children[i];
		if (field == "samples") {
			if (!fv.IsNull()) {
				sel.has_samples = true;
				sel.samples = fv;
			}
		} else if (field == "region") {
			if (!fv.IsNull()) {
				sel.has_region = true;
				sel.region = fv;
			}
		} else if (field == "variants") {
			if (!fv.IsNull()) {
				sel.has_variants = true;
				sel.variants = fv;
			}
		} else if (field == "variant_range") {
			if (!fv.IsNull()) {
				// Half-open [lo, hi) in variant ordinal space; validated at parse time.
				// ResolveSelection wires this into the scanned interval (for ld, the
				// field is rejected up front — see RejectSelectionVariantAxis).
				if (fv.type().id() != LogicalTypeId::LIST) {
					throw InvalidInputException(
					    "selection: variant_range must be a 2-element list [lo, hi) in variant ordinal space (got %s)",
					    fv.type().ToString());
				}
				auto &vr_children = ListValue::GetChildren(fv);
				if (vr_children.size() != 2) {
					throw InvalidInputException(
					    "selection: variant_range must have exactly 2 elements [lo, hi) (got %llu)",
					    static_cast<unsigned long long>(vr_children.size()));
				}
				if (vr_children[0].IsNull() || vr_children[1].IsNull()) {
					throw InvalidInputException("selection: variant_range bounds must not be NULL");
				}
				if (!vr_children[0].type().IsIntegral() || !vr_children[1].type().IsIntegral()) {
					throw InvalidInputException("selection: variant_range bounds must be integers");
				}
				int64_t lo = vr_children[0].GetValue<int64_t>();
				int64_t hi = vr_children[1].GetValue<int64_t>();
				if (lo < 0 || hi < lo) {
					throw InvalidInputException(
					    "selection: variant_range must be a half-open [lo, hi) with 0 <= lo <= hi (got [%lld, %lld))",
					    static_cast<long long>(lo), static_cast<long long>(hi));
				}
				sel.has_variant_range = true;
				sel.variant_range_lo = lo;
				sel.variant_range_hi = hi;
			}
		} else {
			throw InvalidInputException(
			    "selection: unknown field '%s' (expected samples, variants, region, or variant_range)",
			    CompatNameStr(field));
		}
	}
}

Selection ParseSelection(const NamedParamMap &params) {
	Selection sel;

	// --- Lowest precedence: a shared `selection := <STRUCT>` seeds the axes ---
	// A NULL selection means "no selection". Fields are merged below; individual
	// named params then override per-field (most-specific-wins).
	auto selection_it = params.find("selection");
	if (selection_it != params.end() && !selection_it->second.IsNull()) {
		MergeSelectionStruct(selection_it->second, sel);
	}

	// --- Highest precedence: individual named params override per-field ---
	// Unconditional overwrite (not guarded on has_*) so that a present individual
	// param always wins, including an explicit `samples := NULL` (which keeps
	// presence true with a NULL value, so ResolveSampleIndices throws "samples list
	// must not be empty" rather than falling back to the selection's samples).
	//
	// Sample axis: capture the raw Value (ANY); presence is tracked separately from
	// null-ness. Variant axis: keep the region SYMBOLIC (raw Value) — do not parse
	// or resolve it here, so string-format validation and ordinal resolution stay in
	// ResolveSelection at the same logical point the inline code ran them.
	auto samples_it = params.find("samples");
	if (samples_it != params.end()) {
		sel.has_samples = true;
		sel.samples = samples_it->second;
	}

	auto region_it = params.find("region");
	if (region_it != params.end()) {
		sel.has_region = true;
		sel.region = region_it->second;
	}

	// `variants` / `variant_range` have no individual-param surface on the routed
	// functions; they arrive only via the `selection` struct. ResolveSelection then
	// rejects `variants` (arbitrary set, unsupported by a contiguous scan) and wires
	// `variant_range` into the scanned interval.

	return sel;
}

ResolvedSelection ResolveSelection(const Selection &selection, const VariantMetadataIndex &variants,
                                   const SampleInfo *sample_info, uint32_t raw_sample_ct, uint32_t raw_variant_ct,
                                   const string &func_name, bool allow_variant_set) {
	ResolvedSelection resolved;
	resolved.effective_sample_ct = raw_sample_ct;

	// --- sample axis (resolved BEFORE the variant axis to preserve error ordering) ---
	if (selection.has_samples) {
		auto indices = ResolveSampleIndices(selection.samples, raw_sample_ct, sample_info, func_name);

		// Sorted-unique original indices matching pgenlib PgrGet subset ordering,
		// used downstream to align a sex array to effective sample order.
		resolved.subset_sorted = indices;
		std::sort(resolved.subset_sorted.begin(), resolved.subset_sorted.end());
		resolved.subset_sorted.erase(std::unique(resolved.subset_sorted.begin(), resolved.subset_sorted.end()),
		                             resolved.subset_sorted.end());
		resolved.has_subset = true;

		resolved.sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(raw_sample_ct, indices));
		resolved.has_sample_subset = true;
		resolved.effective_sample_ct = resolved.sample_subset->subset_sample_ct;
	}

	// --- variant axis ---
	// These callers scan a single contiguous [start, end) interval, so they support
	// only the contiguous forms of the variant axis. An arbitrary variant SET
	// (`variants` — Tier C) cannot map onto one interval: reject it rather than
	// silently ignore (no false affordance).
	//
	// `allow_variant_set` is the opt-in for a caller that scans a materialized ordinal
	// list and CAN express an arbitrary set (only read_pfile). When it is true, the
	// `variants` field is resolved via the shared ResolveVariantsParameter into a
	// sorted-unique ordinal set, centralizing the arbitrary-SET resolution that used to
	// live inline in read_pfile's bind. `variant_range` is left inactive in this case —
	// read_pfile composes region/variant_range with the set itself (see the law block in
	// plink_common.hpp). The DEFAULT-FALSE throw below is load-bearing: it preserves the
	// analysis functions' rejection + exact error string unchanged.
	if (selection.has_variants) {
		if (!allow_variant_set) {
			throw InvalidInputException("%s: arbitrary variant selection (`variants`) is not supported here; "
			                            "use `region` or `variant_range` for a contiguous range",
			                            func_name);
		}
		resolved.variant_set_sorted = ResolveVariantsParameter(selection.variants, variants, raw_variant_ct, func_name);
		std::sort(resolved.variant_set_sorted.begin(), resolved.variant_set_sorted.end());
		// ResolveVariantsParameter already rejects duplicates, so the sort alone yields a
		// sorted-unique list (no std::unique needed); this matches read_pfile's former
		// per-source variant_set membership exactly.
		resolved.has_variant_set = true;
		return resolved; // variant_range stays inactive; caller owns region/variant_range composition
	}
	if (selection.has_variant_range) {
		// `variant_range := [lo, hi)` is already the contiguous ordinal interval the
		// scan consumes, so wire it straight into VariantRange, clamped to the
		// available variant count (mirrors read_pfile's variant_range handling).
		if (selection.has_region) {
			throw InvalidInputException("%s: specify either `region` or `variant_range`, not both", func_name);
		}
		uint32_t start = static_cast<uint32_t>(std::min<int64_t>(selection.variant_range_lo, raw_variant_ct));
		uint32_t end = static_cast<uint32_t>(std::min<int64_t>(selection.variant_range_hi, raw_variant_ct));
		resolved.variant_range.has_filter = true;
		resolved.variant_range.start_idx = start;
		resolved.variant_range.end_idx = end;
	} else if (selection.has_region) {
		resolved.variant_range = ParseRegion(selection.region.GetValue<string>(), variants, func_name);
	}

	return resolved;
}

void RejectSelectionVariantAxis(const NamedParamMap &params, const string &func_name) {
	// For a function that does NOT take the variant axis from the `selection` struct
	// (currently only plink_ld, whose region grammar differs and which does not
	// support arbitrary/ordinal variant selection): a variant-axis FIELD inside the
	// struct would otherwise be silently dropped. Reject it with a clear error. The
	// individual `region :=` param is untouched — this inspects only the struct.
	auto it = params.find("selection");
	if (it == params.end() || it->second.IsNull()) {
		return;
	}
	const Value &sv = it->second;
	if (sv.type().id() != LogicalTypeId::STRUCT) {
		return; // MergeSelectionStruct (via ParseSelection) reports the type error
	}
	auto &child_types = StructType::GetChildTypes(sv.type());
	auto &children = StructValue::GetChildren(sv);
	for (idx_t i = 0; i < child_types.size(); i++) {
		const auto &field = child_types[i].first; // Identifier on v2.0, string on v1.5; compared below
		if (children[i].IsNull()) {
			continue;
		}
		if (field == "region") {
			throw InvalidInputException("%s: pass `region :=` directly; its region grammar differs and is not "
			                            "taken from the `selection` struct",
			                            func_name);
		}
		if (field == "variants" || field == "variant_range") {
			throw InvalidInputException("%s: `%s` in the `selection` struct is not supported here; "
			                            "use `region :=` directly",
			                            func_name, field);
		}
	}
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

// ---------------------------------------------------------------------------
// Ploidy- and sex-aware statistics (chrX/Y/MT)
// ---------------------------------------------------------------------------

//! Normalize a chromosome label for comparison: lowercase, strip a leading "chr".
static string NormalizeChrom(const string &chrom) {
	string lower = StringUtil::Lower(chrom);
	if (lower.rfind("chr", 0) == 0) {
		lower = lower.substr(3);
	}
	return lower;
}

ParBounds ResolveParBounds(const string &build, const string &func_name) {
	string b = StringUtil::Lower(build);
	// Strip common separators/spaces so "GRCh38", "grch-38", "b38" all match.
	string norm;
	for (char c : b) {
		if (c != '-' && c != '_' && c != ' ' && c != '.') {
			norm.push_back(c);
		}
	}
	ParBounds pb;
	if (norm.empty() || norm == "none") {
		pb.active = false; // rely solely on explicit PAR chromosome codes
		return pb;
	}
	if (norm == "grch38" || norm == "hg38" || norm == "b38" || norm == "38") {
		// GRCh38 PAR (1-based, inclusive): PAR1 X:10001-2781479, PAR2 X:155701383-156030895
		pb.par1_end = 2781479;
		pb.par2_start = 155701383;
		pb.par2_end = 156030895;
		pb.active = true;
		return pb;
	}
	if (norm == "grch37" || norm == "hg19" || norm == "b37" || norm == "37") {
		// GRCh37 PAR (1-based, inclusive): PAR1 X:60001-2699520, PAR2 X:154931044-155260560
		pb.par1_end = 2699520;
		pb.par2_start = 154931044;
		pb.par2_end = 155260560;
		pb.active = true;
		return pb;
	}
	throw InvalidInputException("%s: unrecognized build '%s' (expected 'GRCh38'/'hg38', 'GRCh37'/'hg19', or 'none')",
	                            func_name, build);
}

ChromPloidy ClassifyChromPloidy(const string &chrom, int32_t pos, const ParBounds &par) {
	string c = NormalizeChrom(chrom);
	// Explicit pseudo-autosomal chromosome codes → diploid in both sexes.
	if (c == "par1" || c == "par2" || c == "xy" || c == "25") {
		return ChromPloidy::AUTOSOMAL;
	}
	if (c == "y" || c == "24") {
		return ChromPloidy::CHR_Y;
	}
	if (c == "mt" || c == "m" || c == "26" || c == "chrm") {
		return ChromPloidy::CHR_MT;
	}
	if (c == "x" || c == "23") {
		if (par.active && ((pos > 0 && pos <= par.par1_end) || (pos >= par.par2_start && pos <= par.par2_end))) {
			return ChromPloidy::AUTOSOMAL; // pseudo-autosomal region
		}
		return ChromPloidy::CHR_X;
	}
	return ChromPloidy::AUTOSOMAL;
}

vector<uint8_t> BuildAlignedSex(const SampleInfo &sample_info, const vector<uint32_t> *subset_sorted) {
	if (sample_info.sexes.empty()) {
		return {};
	}
	if (subset_sorted == nullptr) {
		return sample_info.sexes; // full cohort, file order == effective order
	}
	vector<uint8_t> aligned;
	aligned.reserve(subset_sorted->size());
	for (uint32_t idx : *subset_sorted) {
		aligned.push_back(idx < sample_info.sexes.size() ? sample_info.sexes[idx] : uint8_t {0});
	}
	return aligned;
}

SexAwareCounts ComputeSexAwareCounts(const int8_t *geno_bytes, uint32_t sample_ct, ChromPloidy ploidy,
                                     const uint8_t *sex, bool have_sex) {
	SexAwareCounts r;

	// chrX/chrY require sex to stratify ploidy; without it we cannot produce a
	// correct number, so signal NULL rather than silently assuming diploid.
	if ((ploidy == ChromPloidy::CHR_X || ploidy == ChromPloidy::CHR_Y) && !have_sex) {
		r.sex_unavailable = true;
		return r;
	}

	for (uint32_t i = 0; i < sample_ct; i++) {
		int8_t g = geno_bytes[i];
		uint8_t s = have_sex ? sex[i] : uint8_t {0};

		switch (ploidy) {
		case ChromPloidy::CHR_MT: {
			// Haploid for everyone; a heterozygous hardcall is invalid → missing.
			if (g == 0) {
				r.obs_allele_ct += 1;
				r.geno_hom_ref += 1;
			} else if (g == 2) {
				r.obs_allele_ct += 1;
				r.alt_allele_ct += 1;
				r.geno_hom_alt += 1;
			} else {
				r.geno_missing += 1; // g == -9 (missing) or g == 1 (invalid het)
			}
			break;
		}
		case ChromPloidy::CHR_Y: {
			// Males haploid; females and unknown-sex excluded (no Y).
			if (s != 1) {
				r.geno_missing += 1;
				break;
			}
			if (g == 0) {
				r.obs_allele_ct += 1;
				r.geno_hom_ref += 1;
			} else if (g == 2) {
				r.obs_allele_ct += 1;
				r.alt_allele_ct += 1;
				r.geno_hom_alt += 1;
			} else {
				r.geno_missing += 1; // missing or invalid het on haploid male
			}
			break;
		}
		case ChromPloidy::CHR_X: {
			if (s == 2) {
				// Female: diploid. HWE is defined on this stratum.
				if (g == -9) {
					r.geno_missing += 1;
					break;
				}
				r.obs_allele_ct += 2;
				r.alt_allele_ct += static_cast<uint32_t>(g);
				if (g == 0) {
					r.hwe_hom_ref += 1;
					r.geno_hom_ref += 1;
				} else if (g == 1) {
					r.hwe_het += 1;
					r.geno_het += 1;
				} else { // g == 2
					r.hwe_hom_alt += 1;
					r.geno_hom_alt += 1;
				}
			} else if (s == 1) {
				// Male: haploid. Contributes to allele/genotype tallies but NOT to
				// the female-only HWE stratum. Invalid het → missing.
				if (g == 0) {
					r.obs_allele_ct += 1;
					r.geno_hom_ref += 1;
				} else if (g == 2) {
					r.obs_allele_ct += 1;
					r.alt_allele_ct += 1;
					r.geno_hom_alt += 1;
				} else {
					r.geno_missing += 1;
				}
			} else {
				// Unknown sex: cannot determine ploidy → exclude from stratum.
				r.geno_missing += 1;
			}
			break;
		}
		case ChromPloidy::AUTOSOMAL:
		default:
			// Not expected here (callers use the fast diploid path); treat diploid.
			if (g == -9) {
				r.geno_missing += 1;
			} else {
				r.obs_allele_ct += 2;
				r.alt_allele_ct += static_cast<uint32_t>(g);
				if (g == 0) {
					r.hwe_hom_ref += 1;
					r.geno_hom_ref += 1;
				} else if (g == 1) {
					r.hwe_het += 1;
					r.geno_het += 1;
				} else {
					r.hwe_hom_alt += 1;
					r.geno_hom_alt += 1;
				}
			}
			break;
		}
	}

	// HWE is only defined on the diploid stratum of chrX non-PAR (females).
	r.hwe_defined = (ploidy == ChromPloidy::CHR_X || ploidy == ChromPloidy::AUTOSOMAL);
	return r;
}

// ---------------------------------------------------------------------------
// Shared genotype output vector filling
// ---------------------------------------------------------------------------

void FillGenotypeVector(Vector &vec, idx_t row_idx, GenotypeMode mode, uint32_t output_sample_ct,
                        const int8_t *genotype_bytes, const int8_t *phased_pairs, bool include_phased) {
	if (include_phased) {
		if (mode == GenotypeMode::ARRAY) {
			auto &pair_vec = ArrayVector::GetEntry(vec);
			auto &allele_vec = ArrayVector::GetEntry(pair_vec);
			auto *allele_data = CompatFlatDataMutable<int8_t>(allele_vec);
			auto &pair_validity = CompatFlatValidityMutable(pair_vec);

			idx_t pair_base = row_idx * static_cast<idx_t>(output_sample_ct);
			for (idx_t s = 0; s < output_sample_ct; s++) {
				idx_t pair_idx = pair_base + s;
				idx_t allele_base = pair_idx * 2;
				int8_t a1 = phased_pairs[s * 2];
				int8_t a2 = phased_pairs[s * 2 + 1];
				if (a1 == -9) {
					pair_validity.SetInvalid(pair_idx);
					allele_data[allele_base] = 0;
					allele_data[allele_base + 1] = 0;
				} else {
					allele_data[allele_base] = a1;
					allele_data[allele_base + 1] = a2;
				}
			}
		} else {
			// LIST mode
			auto list_offset = ListVector::GetListSize(vec);
			ListVector::Reserve(vec, list_offset + output_sample_ct);
			auto &pair_vec = ListVector::GetEntry(vec);
			auto &allele_vec = ArrayVector::GetEntry(pair_vec);
			auto *allele_data = CompatFlatDataMutable<int8_t>(allele_vec);
			auto &pair_validity = CompatFlatValidityMutable(pair_vec);

			for (idx_t s = 0; s < output_sample_ct; s++) {
				idx_t pair_idx = list_offset + s;
				idx_t allele_base = pair_idx * 2;
				int8_t a1 = phased_pairs[s * 2];
				int8_t a2 = phased_pairs[s * 2 + 1];
				if (a1 == -9) {
					pair_validity.SetInvalid(pair_idx);
					allele_data[allele_base] = 0;
					allele_data[allele_base + 1] = 0;
				} else {
					allele_data[allele_base] = a1;
					allele_data[allele_base + 1] = a2;
				}
			}

			auto *list_data = CompatFlatDataMutable<list_entry_t>(vec);
			list_data[row_idx].offset = list_offset;
			list_data[row_idx].length = output_sample_ct;
			ListVector::SetListSize(vec, list_offset + output_sample_ct);
		}
	} else {
		// Unphased
		if (mode == GenotypeMode::ARRAY) {
			auto array_size = static_cast<idx_t>(output_sample_ct);
			auto &child = ArrayVector::GetEntry(vec);
			auto *child_data = CompatFlatDataMutable<int8_t>(child);
			auto &child_validity = CompatFlatValidityMutable(child);

			idx_t base = row_idx * array_size;
			for (idx_t s = 0; s < array_size; s++) {
				int8_t geno = genotype_bytes[s];
				if (geno == -9) {
					child_validity.SetInvalid(base + s);
					child_data[base + s] = 0;
				} else {
					child_data[base + s] = geno;
				}
			}
		} else {
			// LIST mode
			auto list_offset = ListVector::GetListSize(vec);
			ListVector::Reserve(vec, list_offset + output_sample_ct);
			auto &child = ListVector::GetEntry(vec);
			auto *child_data = CompatFlatDataMutable<int8_t>(child);
			auto &child_validity = CompatFlatValidityMutable(child);
			for (idx_t s = 0; s < output_sample_ct; s++) {
				int8_t geno = genotype_bytes[s];
				if (geno == -9) {
					child_validity.SetInvalid(list_offset + s);
					child_data[list_offset + s] = 0;
				} else {
					child_data[list_offset + s] = geno;
				}
			}
			auto *list_data = CompatFlatDataMutable<list_entry_t>(vec);
			list_data[row_idx].offset = list_offset;
			list_data[row_idx].length = output_sample_ct;
			ListVector::SetListSize(vec, list_offset + output_sample_ct);
		}
	}
}

// ---------------------------------------------------------------------------
// pfile catalog grid
// ---------------------------------------------------------------------------

//! A v1 catalog as a grid: concat_dim "variants" is N variant blocks x 1 sample block, "samples"
//! is 1 x K, and a single-fileset catalog is 1 x 1. The parser has already folded the shared
//! companion (.psam or .pvar) into every entry and resolved every path.
static void LoadV1CatalogGrid(FileSystem &fs, const string &content, CatalogGrid &grid, const string &func_name) {
	auto cat = ParseV1Catalog(fs, grid.path, content);
	auto files = cat.files;
	std::stable_sort(files.begin(), files.end(),
	                 [](const ShardEntry &a, const ShardEntry &b) { return a.partition_start < b.partition_start; });
	const bool by_samples = cat.concat_dim == "samples";
	grid.has_variant_extent = !cat.is_single && !by_samples;
	grid.has_sample_extent = !cat.is_single && by_samples;
	grid.sample_block_ct = by_samples ? files.size() : 1;
	for (idx_t i = 0; i < files.size(); i++) {
		MatrixCell cell;
		cell.prefix = files[i].prefix;
		cell.pvar = files[i].pvar;
		cell.psam = files[i].psam;
		cell.variant_block = by_samples ? 0 : static_cast<int64_t>(i);
		cell.sample_block = by_samples ? static_cast<int64_t>(i) : 0;
		auto &start = by_samples ? cell.sample_start : cell.variant_start;
		auto &stop = by_samples ? cell.sample_stop : cell.variant_stop;
		start = files[i].partition_start;
		stop = files[i].partition_stop;
		grid.cells.push_back(std::move(cell));
	}
	if (grid.cells.empty()) {
		throw InvalidInputException("%s: catalog '%s' lists no files", func_name, grid.path);
	}
	if (by_samples || cat.is_single) {
		grid.variant_ranges.resize(1); // one variant block, no declared range: never pruned
	} else {
		grid.variant_ranges = files; // each shard's own `range` (CHROM list, POS extent)
	}
}

CatalogGrid LoadCatalogGrid(ClientContext &context, FileSystem &fs, const string &path, const string &func_name) {
	CatalogGrid grid;
	grid.path = ResolveExistingPath(context, fs, path);
	if (grid.path.empty()) {
		throw InvalidInputException("%s: catalog file not found: '%s'", func_name, path);
	}
	string content;
	{
		auto handle = fs.OpenFile(grid.path, FileFlags::FILE_FLAGS_READ);
		content.resize(handle->GetFileSize());
		if (!content.empty()) {
			handle->Read(reinterpret_cast<void *>(&content[0]), content.size());
		}
	}
	if (!IsMatrixCatalog(content)) {
		LoadV1CatalogGrid(fs, content, grid, func_name);
		grid.declared_variant_block_ct = grid.variant_ranges.size();
		grid.declared_sample_block_ct = grid.sample_block_ct;
		return grid;
	}
	auto matrix = ParseMatrixCatalog(fs, grid.path, content);
	grid.cells = BuildCells(matrix, fs, /*include_holes=*/true);
	for (auto &cell : grid.cells) {
		if (!cell.present) {
			string where;
			for (auto &kv : cell.keys) {
				where += (where.empty() ? "" : ", ") + kv.first + "=" + kv.second;
			}
			throw InvalidInputException("%s: catalog '%s' has a hole (an absent cell) at {%s}. A 2-D read of "
			                            "a matrix with holes is not supported yet.",
			                            func_name, path, where);
		}
	}
	std::stable_sort(grid.cells.begin(), grid.cells.end(), [](const MatrixCell &a, const MatrixCell &b) {
		return a.variant_block != b.variant_block ? a.variant_block < b.variant_block : a.sample_block < b.sample_block;
	});
	grid.sample_block_ct = matrix.samples.blocks.size();
	const idx_t variant_block_ct = matrix.variants.blocks.size();
	if (grid.cells.empty() || grid.sample_block_ct == 0 ||
	    grid.cells.size() != variant_block_ct * grid.sample_block_ct) {
		throw InvalidInputException("%s: catalog '%s' does not describe a full grid (%llu cells for %llu "
		                            "variant blocks x %llu sample blocks)",
		                            func_name, path, static_cast<unsigned long long>(grid.cells.size()),
		                            static_cast<unsigned long long>(variant_block_ct),
		                            static_cast<unsigned long long>(grid.sample_block_ct));
	}
	for (idx_t i = 0; i < grid.cells.size(); i++) {
		auto &cell = grid.cells[i];
		if (static_cast<idx_t>(cell.variant_block) != i / grid.sample_block_ct ||
		    static_cast<idx_t>(cell.sample_block) != i % grid.sample_block_ct) {
			throw InvalidInputException("%s: catalog '%s' does not describe a full grid", func_name, path);
		}
	}
	// Each variant block's declared range: its CHROM enum value, and the inclusive POS bounds of an
	// ordinal_interval key over POS. A block that declares neither has no range and is never pruned.
	const AxisKey *contig_key = matrix.variants.ContigKey();
	const AxisKey *pos_key = nullptr;
	for (auto &k : matrix.variants.keys) {
		if (k.type == AxisKeyType::ORDINAL_INTERVAL && k.over == "POS") {
			pos_key = &k;
			break;
		}
	}
	for (auto &block : matrix.variants.blocks) {
		ShardEntry r;
		if (contig_key) {
			r.chrom.push_back(block.ValueOf(contig_key->name));
			r.has_range = true;
		}
		const BlockBound *b = pos_key ? block.FindBound(pos_key->name) : nullptr;
		if (b) {
			r.has_pos = true;
			r.pos_min = b->lo;
			r.pos_max = b->hi;
			r.has_range = true;
		}
		grid.variant_ranges.push_back(std::move(r));
	}
	grid.declared_variant_block_ct = variant_block_ct;
	grid.declared_sample_block_ct = grid.sample_block_ct;
	return grid;
}

} // namespace duckdb
