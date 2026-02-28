#include "plink_common.hpp"

namespace duckdb {

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
	                        static_cast<unsigned long long>(field_idx),
	                        static_cast<unsigned long long>(vidx));
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
		throw InternalException("VariantMetadataIndex::GetPos: invalid POS value '%s' for variant %llu",
		                        field.c_str(), static_cast<unsigned long long>(vidx));
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

} // namespace duckdb
