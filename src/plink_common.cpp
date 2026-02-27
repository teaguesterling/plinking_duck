#include "plink_common.hpp"

namespace duckdb {

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
// Metadata loading
// ---------------------------------------------------------------------------

VariantMetadata LoadVariantMetadata(ClientContext &context, const string &path, const string &func_name) {
	// Read the file once into memory, then parse header and data from the
	// in-memory lines. Avoids the double-read that would occur if we called
	// ParsePvarHeader (which opens the file separately) before ReadFileLines.
	auto lines = ReadFileLines(context, path);

	if (lines.empty()) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is empty", func_name, path);
	}

	// Detect format and parse header from in-memory lines
	idx_t skip_lines = 0;
	bool is_bim = false;
	vector<string> column_names;

	// Skip ## comment/meta lines
	while (skip_lines < lines.size()) {
		auto &line = lines[skip_lines];
		if (line.empty() || (line.size() >= 2 && line[0] == '#' && line[1] == '#')) {
			skip_lines++;
			continue;
		}
		break;
	}

	if (skip_lines >= lines.size()) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' contains no header or data", func_name, path);
	}

	auto &header_line = lines[skip_lines];
	if (header_line.size() >= 6 && header_line.substr(0, 6) == "#CHROM") {
		// .pvar format: parse column names from the #CHROM header line
		is_bim = false;
		auto fields = SplitTabLine(header_line.substr(1)); // strip leading '#'
		column_names = std::move(fields);
		skip_lines++; // skip the header line
	} else {
		// Legacy .bim format: fixed 6-column schema, no header to skip
		// Normalized order: CHROM, POS, ID, REF, ALT, CM
		is_bim = true;
		column_names = {"CHROM", "POS", "ID", "REF", "ALT", "CM"};
	}

	// Find column indices for the 5 core fields
	idx_t chrom_idx = DConstants::INVALID_INDEX;
	idx_t pos_idx = DConstants::INVALID_INDEX;
	idx_t id_idx = DConstants::INVALID_INDEX;
	idx_t ref_idx = DConstants::INVALID_INDEX;
	idx_t alt_idx = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < column_names.size(); i++) {
		const auto &name = column_names[i];
		if (name == "CHROM") {
			chrom_idx = i;
		} else if (name == "POS") {
			pos_idx = i;
		} else if (name == "ID") {
			id_idx = i;
		} else if (name == "REF") {
			ref_idx = i;
		} else if (name == "ALT") {
			alt_idx = i;
		}
	}

	if (chrom_idx == DConstants::INVALID_INDEX || pos_idx == DConstants::INVALID_INDEX ||
	    id_idx == DConstants::INVALID_INDEX || ref_idx == DConstants::INVALID_INDEX ||
	    alt_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("%s: .pvar/.bim file '%s' is missing required columns "
		                            "(need CHROM, POS, ID, REF, ALT)",
		                            func_name, path);
	}

	VariantMetadata meta;
	for (idx_t line_idx = skip_lines; line_idx < lines.size(); line_idx++) {
		auto &line = lines[line_idx];
		if (line.empty()) {
			continue;
		}

		auto fields = is_bim ? SplitWhitespaceLine(line) : SplitTabLine(line);

		// .bim files: CHROM(0) ID(1) CM(2) POS(3) ALT(4) REF(5)
		// Normalize to: CHROM(0) POS(1) ID(2) REF(3) ALT(4) CM(5)
		vector<string> *source = &fields;
		vector<string> normalized;
		if (is_bim) {
			if (fields.size() < 6) {
				throw InvalidInputException("%s: .bim file '%s' has line with %llu fields, expected 6", func_name, path,
				                            static_cast<unsigned long long>(fields.size()));
			}
			normalized = {fields[0], fields[3], fields[1], fields[5], fields[4], fields[2]};
			source = &normalized;
		}

		auto &src = *source;
		if (chrom_idx >= src.size() || pos_idx >= src.size() || id_idx >= src.size() || ref_idx >= src.size() ||
		    alt_idx >= src.size()) {
			throw InvalidInputException("%s: .pvar/.bim file '%s' has line with too few fields", func_name, path);
		}

		meta.chroms.push_back(src[chrom_idx]);

		char *end;
		errno = 0;
		long pos_val = std::strtol(src[pos_idx].c_str(), &end, 10);
		if (end == src[pos_idx].c_str() || *end != '\0' || errno != 0) {
			throw InvalidInputException("%s: invalid POS value '%s' in '%s'", func_name, src[pos_idx], path);
		}
		meta.positions.push_back(static_cast<int32_t>(pos_val));

		meta.ids.push_back(src[id_idx] == "." ? "" : src[id_idx]);
		meta.refs.push_back(src[ref_idx]);
		meta.alts.push_back(src[alt_idx]);
	}

	meta.variant_ct = meta.chroms.size();
	return meta;
}

// ---------------------------------------------------------------------------
// Sample parameter resolution
// ---------------------------------------------------------------------------

vector<uint32_t> ResolveSampleIndices(const Value &samples_val, uint32_t raw_sample_ct,
                                      const SampleInfo *sample_info, const string &func_name) {
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

VariantRange ParseRegion(const string &region_str, const VariantMetadata &variants, const string &func_name) {
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

	// Scan variant metadata for matching range
	// Assumes sorted by (CHROM, POS) per PLINK format
	VariantRange range;
	range.has_filter = true;

	bool found_start = false;
	for (uint32_t i = 0; i < static_cast<uint32_t>(variants.variant_ct); i++) {
		if (variants.chroms[i] == chrom && variants.positions[i] >= static_cast<int32_t>(start_pos) &&
		    variants.positions[i] <= static_cast<int32_t>(end_pos)) {
			if (!found_start) {
				range.start_idx = i;
				found_start = true;
			}
			range.end_idx = i + 1;
		}
	}

	return range;
}

} // namespace duckdb
