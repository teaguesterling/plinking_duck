#include "pfile_reader.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"

#include "duckdb/common/file_system.hpp"

#include <pgenlib_read.h>
#include <pgenlib_ffi_support.h>

#include <algorithm>
#include <atomic>
#include <cstring>
#include <unordered_set>

namespace duckdb {

// ---------------------------------------------------------------------------
// RAII wrapper for cache-aligned allocations (same pattern as pgen_reader)
// ---------------------------------------------------------------------------

struct PfileAlignedBuffer {
	void *ptr = nullptr;

	~PfileAlignedBuffer() {
		if (ptr) {
			plink2::aligned_free(ptr);
		}
	}

	PfileAlignedBuffer() = default;
	PfileAlignedBuffer(const PfileAlignedBuffer &) = delete;
	PfileAlignedBuffer &operator=(const PfileAlignedBuffer &) = delete;
	PfileAlignedBuffer(PfileAlignedBuffer &&other) noexcept : ptr(other.ptr) {
		other.ptr = nullptr;
	}
	PfileAlignedBuffer &operator=(PfileAlignedBuffer &&other) noexcept {
		if (this != &other) {
			if (ptr) {
				plink2::aligned_free(ptr);
			}
			ptr = other.ptr;
			other.ptr = nullptr;
		}
		return *this;
	}

	void Allocate(uintptr_t size) {
		if (plink2::cachealigned_malloc(size, &ptr)) {
			throw IOException("read_pfile: failed to allocate %llu bytes of aligned memory",
			                  static_cast<unsigned long long>(size));
		}
	}

	template <typename T>
	T *As() {
		return static_cast<T *>(ptr);
	}
};

// ---------------------------------------------------------------------------
// Variant metadata (same as pgen_reader, duplicated to avoid coupling)
// ---------------------------------------------------------------------------

struct PfileVariantMetadata {
	vector<string> chroms;
	vector<int32_t> positions;
	vector<string> ids;
	vector<string> refs;
	vector<string> alts;
	idx_t variant_ct = 0;

	//! Build ID→index map for variant filtering by name
	unordered_map<string, uint32_t> id_to_idx;

	void BuildIdMap() {
		for (idx_t i = 0; i < variant_ct; i++) {
			if (!ids[i].empty()) {
				id_to_idx[ids[i]] = static_cast<uint32_t>(i);
			}
		}
	}
};

// ---------------------------------------------------------------------------
// File reading and variant metadata loading (reused from pgen_reader pattern)
// ---------------------------------------------------------------------------

static vector<string> PfileReadFileLines(ClientContext &context, const string &path) {
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

static vector<string> PfileSplitTabLine(const string &line) {
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

static vector<string> PfileSplitWhitespaceLine(const string &line) {
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

static PfileVariantMetadata LoadPfileVariantMetadata(ClientContext &context, const string &path) {
	auto header_info = ParsePvarHeader(context, path);

	idx_t chrom_idx = DConstants::INVALID_INDEX;
	idx_t pos_idx = DConstants::INVALID_INDEX;
	idx_t id_idx = DConstants::INVALID_INDEX;
	idx_t ref_idx = DConstants::INVALID_INDEX;
	idx_t alt_idx = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < header_info.column_names.size(); i++) {
		const auto &name = header_info.column_names[i];
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
		throw InvalidInputException("read_pfile: .pvar/.bim file '%s' is missing required columns "
		                            "(need CHROM, POS, ID, REF, ALT)",
		                            path);
	}

	auto lines = PfileReadFileLines(context, path);

	PfileVariantMetadata meta;
	for (idx_t line_idx = header_info.skip_lines; line_idx < lines.size(); line_idx++) {
		auto &line = lines[line_idx];
		if (line.empty()) {
			continue;
		}

		auto fields = header_info.is_bim ? PfileSplitWhitespaceLine(line) : PfileSplitTabLine(line);

		vector<string> *source = &fields;
		vector<string> normalized;
		if (header_info.is_bim) {
			if (fields.size() < 6) {
				throw InvalidInputException("read_pfile: .bim file '%s' has line with %llu fields, expected 6", path,
				                            static_cast<unsigned long long>(fields.size()));
			}
			normalized = {fields[0], fields[3], fields[1], fields[5], fields[4], fields[2]};
			source = &normalized;
		}

		auto &src = *source;
		if (chrom_idx >= src.size() || pos_idx >= src.size() || id_idx >= src.size() || ref_idx >= src.size() ||
		    alt_idx >= src.size()) {
			throw InvalidInputException("read_pfile: .pvar/.bim file '%s' has line with too few fields", path);
		}

		meta.chroms.push_back(src[chrom_idx]);

		char *end;
		errno = 0;
		long pos_val = std::strtol(src[pos_idx].c_str(), &end, 10);
		if (end == src[pos_idx].c_str() || *end != '\0' || errno != 0) {
			throw InvalidInputException("read_pfile: invalid POS value '%s' in '%s'", src[pos_idx], path);
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
// Companion file discovery
// ---------------------------------------------------------------------------

static string PfileReplaceExtension(const string &path, const string &new_ext) {
	auto dot = path.rfind('.');
	if (dot == string::npos) {
		return path + new_ext;
	}
	return path.substr(0, dot) + new_ext;
}

static string PfileFindCompanionFile(FileSystem &fs, const string &base_path, const vector<string> &extensions) {
	for (auto &ext : extensions) {
		auto candidate = PfileReplaceExtension(base_path, ext);
		if (fs.FileExists(candidate)) {
			return candidate;
		}
	}
	return "";
}

// ---------------------------------------------------------------------------
// Region parsing
// ---------------------------------------------------------------------------

struct RegionFilter {
	string chrom;
	int64_t start = 0;       // 1-based, inclusive
	int64_t end = INT64_MAX; // 1-based, inclusive
	bool active = false;
};

//! Parse a region string in the format "chr", "chr:start-end", or "chr:start-"
static RegionFilter ParseRegion(const string &region_str) {
	RegionFilter region;
	region.active = true;

	auto colon_pos = region_str.find(':');
	if (colon_pos == string::npos) {
		// Chromosome-only filter
		region.chrom = region_str;
		return region;
	}

	region.chrom = region_str.substr(0, colon_pos);
	if (region.chrom.empty()) {
		throw InvalidInputException("read_pfile: invalid region format '%s' (empty chromosome)", region_str);
	}

	auto range_str = region_str.substr(colon_pos + 1);
	auto dash_pos = range_str.find('-');
	if (dash_pos == string::npos) {
		throw InvalidInputException("read_pfile: invalid region format '%s' (expected chr:start-end)", region_str);
	}

	auto start_str = range_str.substr(0, dash_pos);
	auto end_str = range_str.substr(dash_pos + 1);

	if (start_str.empty()) {
		throw InvalidInputException("read_pfile: invalid region format '%s' (empty start position)", region_str);
	}

	char *parse_end;
	errno = 0;
	region.start = std::strtol(start_str.c_str(), &parse_end, 10);
	if (parse_end == start_str.c_str() || *parse_end != '\0' || errno != 0 || region.start < 0) {
		throw InvalidInputException("read_pfile: invalid region start '%s' in '%s'", start_str, region_str);
	}

	if (!end_str.empty()) {
		errno = 0;
		region.end = std::strtol(end_str.c_str(), &parse_end, 10);
		if (parse_end == end_str.c_str() || *parse_end != '\0' || errno != 0 || region.end < 0) {
			throw InvalidInputException("read_pfile: invalid region end '%s' in '%s'", end_str, region_str);
		}
	}

	if (region.start > region.end) {
		throw InvalidInputException("read_pfile: region start (%lld) > end (%lld) in '%s'",
		                            static_cast<long long>(region.start), static_cast<long long>(region.end),
		                            region_str);
	}

	return region;
}

// ---------------------------------------------------------------------------
// Sample metadata for tidy mode output
// ---------------------------------------------------------------------------

//! Extended sample info for tidy mode: includes all psam columns.
struct PfileSampleMetadata {
	PsamHeaderInfo header;
	//! All data rows from the .psam/.fam, in file order.
	//! Each inner vector has one field per column.
	vector<vector<string>> rows;

	//! Index of SEX column in the header (DConstants::INVALID_INDEX if absent)
	idx_t sex_col_idx = DConstants::INVALID_INDEX;
	//! Indices of PAT/MAT columns for special missing-value handling
	vector<idx_t> parent_col_indices;
};

//! Missing value sentinels (same as psam_reader)
static bool PfileIsMissingValue(const string &val) {
	return val.empty() || val == "." || val == "NA" || val == "na";
}

//! Load full sample metadata including all columns for tidy mode output.
static PfileSampleMetadata LoadPfileSampleMetadata(ClientContext &context, const string &path) {
	auto lines = PfileReadFileLines(context, path);

	PfileSampleMetadata meta;

	// Parse header
	if (lines.empty()) {
		throw IOException("read_pfile: .psam/.fam file '%s' is empty", path);
	}

	meta.header = ParsePsamHeader(context, path);

	// Find special column indices
	for (idx_t i = 0; i < meta.header.column_names.size(); i++) {
		if (meta.header.column_names[i] == "SEX") {
			meta.sex_col_idx = i;
		} else if (meta.header.column_names[i] == "PAT" || meta.header.column_names[i] == "MAT") {
			meta.parent_col_indices.push_back(i);
		}
	}

	// Data lines start at index 1 for .psam (skip header), 0 for .fam
	idx_t data_start = (meta.header.format != PsamFormat::FAM) ? 1 : 0;

	for (idx_t i = data_start; i < lines.size(); i++) {
		auto &line = lines[i];
		if (line.empty()) {
			continue;
		}

		vector<string> fields;
		if (meta.header.format == PsamFormat::FAM) {
			fields = PfileSplitWhitespaceLine(line);
		} else {
			fields = PfileSplitTabLine(line);
		}

		meta.rows.push_back(std::move(fields));
	}

	return meta;
}

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PfileBindData : public TableFunctionData {
	string pgen_path;
	string pvar_path;
	string psam_path;

	// Variant metadata
	PfileVariantMetadata variants;

	// Sample metadata (basic: for IID lookups and default mode)
	SampleInfo sample_info;

	// Full sample metadata (for tidy mode column output)
	PfileSampleMetadata sample_metadata;

	// pgenlib header info
	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Mode
	bool tidy_mode = false;

	// Options
	bool include_dosages = false;
	bool include_phased = false;

	// Sample subsetting
	bool has_sample_subset = false;
	vector<uint32_t> sample_indices; // 0-based indices into .pgen sample order
	uint32_t subset_sample_ct = 0;

	// Region filtering
	RegionFilter region;

	// Variant filtering
	bool has_variant_filter = false;
	vector<uint32_t> variant_indices; // sorted indices of variants to include

	// Effective variant list (after region + variant filter intersection)
	// If empty and no filters active, scan all variants sequentially.
	// If non-empty, these are the specific variant indices to scan.
	bool has_effective_variant_list = false;
	vector<uint32_t> effective_variant_indices;

	//! Number of output samples (after subsetting)
	uint32_t OutputSampleCt() const {
		return has_sample_subset ? subset_sample_ct : static_cast<uint32_t>(sample_info.sample_ct);
	}

	//! Number of effective variants (after filtering)
	uint32_t EffectiveVariantCt() const {
		if (has_effective_variant_list) {
			return static_cast<uint32_t>(effective_variant_indices.size());
		}
		return raw_variant_ct;
	}

	// --- Default mode column indices ---
	static constexpr idx_t CHROM_COL = 0;
	static constexpr idx_t POS_COL = 1;
	static constexpr idx_t ID_COL = 2;
	static constexpr idx_t REF_COL = 3;
	static constexpr idx_t ALT_COL = 4;
	static constexpr idx_t GENOTYPES_COL = 5;

	// --- Tidy mode column layout ---
	// Variant metadata: CHROM(0), POS(1), ID(2), REF(3), ALT(4)
	// Sample metadata: columns from .psam start at index 5
	// Genotype: after all sample metadata columns
	idx_t tidy_sample_col_start = 5;
	idx_t tidy_genotype_col = 0; // computed in bind
	idx_t tidy_total_cols = 0;   // computed in bind

	// Mapping from tidy output sample column index (relative to tidy_sample_col_start)
	// to psam file column index
	vector<idx_t> tidy_sample_col_to_psam_col;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PfileGlobalState : public GlobalTableFunctionState {
	// For default mode: atomic variant counter
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t total_variants = 0;

	// For tidy mode: sequential scanning
	std::atomic<uint32_t> next_effective_variant_pos {0};
	uint32_t total_effective_variants = 0;

	// Projection
	bool need_genotypes = false;
	vector<column_t> column_ids;

	idx_t MaxThreads() const override {
		// Tidy mode must be single-threaded: the state machine tracks
		// current_variant and current_sample, which are not thread-safe.
		// Default mode can parallelize across variants.
		return 1;
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PfileLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	PfileAlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	PfileAlignedBuffer pgr_alloc_buf;
	PfileAlignedBuffer genovec_buf;
	PfileAlignedBuffer sample_include_buf;
	PfileAlignedBuffer cumulative_popcounts_buf;

	vector<int8_t> genotype_bytes;

	plink2::PgrSampleSubsetIndex pssi;

	bool initialized = false;

	// Tidy mode cursor state
	uint32_t tidy_current_variant_pos = 0; // index into effective_variant_indices
	uint32_t tidy_current_sample = 0;      // sample index within current variant
	bool tidy_variant_loaded = false;      // genotypes decoded for current variant
	bool tidy_done = false;

	~PfileLocalState() {
		if (initialized) {
			plink2::PglErr reterr = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr, &reterr);
			plink2::CleanupPgfi(&pgfi, &reterr);
		}
	}
};

// ---------------------------------------------------------------------------
// Bind function
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PfileBind(ClientContext &context, TableFunctionBindInput &input,
                                          vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PfileBindData>();

	// --- Resolve file paths ---
	// First positional argument is the prefix (or empty if only named params)
	string prefix;
	if (!input.inputs.empty()) {
		prefix = input.inputs[0].GetValue<string>();
	}

	auto &fs = FileSystem::GetFileSystem(context);

	// Named parameters can override individual paths
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pgen") {
			bind_data->pgen_path = kv.second.GetValue<string>();
		} else if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "tidy") {
			bind_data->tidy_mode = kv.second.GetValue<bool>();
		} else if (kv.first == "dosages") {
			bind_data->include_dosages = kv.second.GetValue<bool>();
		} else if (kv.first == "phased") {
			bind_data->include_phased = kv.second.GetValue<bool>();
		} else if (kv.first == "region") {
			bind_data->region = ParseRegion(kv.second.GetValue<string>());
		}
		// samples and variants handled after pgenlib init
	}

	if (bind_data->include_dosages) {
		throw NotImplementedException("read_pfile: dosages support is not yet implemented");
	}
	if (bind_data->include_phased) {
		throw NotImplementedException("read_pfile: phased support is not yet implemented");
	}

	// Discover files from prefix if not explicitly provided
	if (bind_data->pgen_path.empty() && !prefix.empty()) {
		// Try prefix.pgen first, then prefix as-is if it already has an extension
		string candidate = prefix + ".pgen";
		if (fs.FileExists(candidate)) {
			bind_data->pgen_path = candidate;
		} else if (fs.FileExists(prefix)) {
			bind_data->pgen_path = prefix;
		} else {
			throw InvalidInputException("read_pfile: cannot find .pgen file for prefix '%s' (tried '%s')", prefix,
			                            candidate);
		}
	}

	if (bind_data->pgen_path.empty()) {
		throw InvalidInputException("read_pfile: no .pgen file path provided");
	}

	if (bind_data->pvar_path.empty()) {
		if (!prefix.empty()) {
			bind_data->pvar_path = PfileFindCompanionFile(fs, prefix + ".pgen", {".pvar", ".bim"});
			if (bind_data->pvar_path.empty()) {
				// Also try prefix-based discovery
				for (auto &ext : {".pvar", ".bim"}) {
					auto candidate = prefix + ext;
					if (fs.FileExists(candidate)) {
						bind_data->pvar_path = candidate;
						break;
					}
				}
			}
		}
		if (bind_data->pvar_path.empty()) {
			bind_data->pvar_path = PfileFindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		}
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("read_pfile: cannot find .pvar or .bim file for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            prefix.empty() ? bind_data->pgen_path : prefix);
		}
	}

	if (bind_data->psam_path.empty()) {
		if (!prefix.empty()) {
			for (auto &ext : {".psam", ".fam"}) {
				auto candidate = prefix + ext;
				if (fs.FileExists(candidate)) {
					bind_data->psam_path = candidate;
					break;
				}
			}
		}
		if (bind_data->psam_path.empty()) {
			bind_data->psam_path = PfileFindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		}
		if (bind_data->psam_path.empty()) {
			throw InvalidInputException("read_pfile: cannot find .psam or .fam file for '%s' "
			                            "(use psam := 'path' to specify explicitly)",
			                            prefix.empty() ? bind_data->pgen_path : prefix);
		}
	}

	// --- Initialize pgenlib (Phase 1) to get counts ---
	plink2::PgenFileInfo pgfi;
	plink2::PreinitPgfi(&pgfi);

	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr, UINT32_MAX, UINT32_MAX,
	                                            &header_ctrl, &pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("read_pfile: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;

	// Phase 2
	PfileAlignedBuffer pgfi_alloc;
	if (pgfi_alloc_cacheline_ct > 0) {
		pgfi_alloc.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	}

	uint32_t max_vrec_width = 0;
	uintptr_t pgr_alloc_cacheline_ct = 0;

	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, pgfi.raw_variant_ct, &max_vrec_width, &pgfi,
	                             pgfi_alloc.As<unsigned char>(), &pgr_alloc_cacheline_ct, errstr_buf);

	plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
	plink2::CleanupPgfi(&pgfi, &cleanup_err);

	if (err != plink2::kPglRetSuccess) {
		throw IOException("read_pfile: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadPfileVariantMetadata(context, bind_data->pvar_path);

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("read_pfile: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info ---
	bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);

	if (bind_data->sample_info.sample_ct != bind_data->raw_sample_ct) {
		throw InvalidInputException("read_pfile: sample count mismatch: .pgen has %u samples, "
		                            ".psam/.fam '%s' has %llu samples",
		                            bind_data->raw_sample_ct, bind_data->psam_path,
		                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
	}

	// Load full sample metadata for tidy mode
	if (bind_data->tidy_mode) {
		bind_data->sample_metadata = LoadPfileSampleMetadata(context, bind_data->psam_path);
	}

	// --- Process samples parameter ---
	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		auto &samples_val = samples_it->second;
		auto &child_type = ListType::GetChildType(samples_val.type());

		auto &children_check = ListValue::GetChildren(samples_val);
		if (children_check.empty()) {
			throw InvalidInputException("read_pfile: samples list must not be empty");
		}

		if (child_type.id() == LogicalTypeId::INTEGER || child_type.id() == LogicalTypeId::BIGINT) {
			auto &children = ListValue::GetChildren(samples_val);
			for (auto &child : children) {
				int64_t idx = child.GetValue<int64_t>();
				if (idx < 0 || static_cast<uint32_t>(idx) >= bind_data->raw_sample_ct) {
					throw InvalidInputException("read_pfile: sample index %lld out of range (sample count: %u)",
					                            static_cast<long long>(idx), bind_data->raw_sample_ct);
				}
				bind_data->sample_indices.push_back(static_cast<uint32_t>(idx));
			}
		} else if (child_type.id() == LogicalTypeId::VARCHAR) {
			auto &children = ListValue::GetChildren(samples_val);
			for (auto &child : children) {
				auto iid = child.GetValue<string>();
				auto it = bind_data->sample_info.iid_to_idx.find(iid);
				if (it == bind_data->sample_info.iid_to_idx.end()) {
					throw InvalidInputException("read_pfile: sample '%s' not found in .psam", iid);
				}
				bind_data->sample_indices.push_back(static_cast<uint32_t>(it->second));
			}
		} else {
			throw InvalidInputException("read_pfile: samples parameter must be LIST(VARCHAR) or LIST(INTEGER)");
		}

		// Validate no duplicates
		{
			std::unordered_set<uint32_t> seen;
			for (auto idx : bind_data->sample_indices) {
				if (!seen.insert(idx).second) {
					throw InvalidInputException("read_pfile: duplicate sample index %u in samples list", idx);
				}
			}
		}

		// Sort sample_indices to match pgenlib output order.
		// PgrGet with sample_include returns genotypes in ascending bit order
		// regardless of the order indices were added. Sorting ensures that
		// genotype_bytes[i] corresponds to sample_indices[i] in tidy mode.
		std::sort(bind_data->sample_indices.begin(), bind_data->sample_indices.end());

		bind_data->has_sample_subset = true;
		bind_data->subset_sample_ct = static_cast<uint32_t>(bind_data->sample_indices.size());
	}

	// --- Process variants parameter ---
	auto variants_it = input.named_parameters.find("variants");
	if (variants_it != input.named_parameters.end()) {
		auto &variants_val = variants_it->second;
		auto &child_type = ListType::GetChildType(variants_val.type());

		auto &children_check = ListValue::GetChildren(variants_val);
		if (children_check.empty()) {
			throw InvalidInputException("read_pfile: variants list must not be empty");
		}

		if (child_type.id() == LogicalTypeId::INTEGER || child_type.id() == LogicalTypeId::BIGINT) {
			auto &children = ListValue::GetChildren(variants_val);
			for (auto &child : children) {
				int64_t idx = child.GetValue<int64_t>();
				if (idx < 0 || static_cast<uint32_t>(idx) >= bind_data->raw_variant_ct) {
					throw InvalidInputException("read_pfile: variant index %lld out of range (variant count: %u)",
					                            static_cast<long long>(idx), bind_data->raw_variant_ct);
				}
				bind_data->variant_indices.push_back(static_cast<uint32_t>(idx));
			}
		} else if (child_type.id() == LogicalTypeId::VARCHAR) {
			// Build ID map if not already built
			bind_data->variants.BuildIdMap();

			auto &children = ListValue::GetChildren(variants_val);
			for (auto &child : children) {
				auto vid = child.GetValue<string>();
				auto it = bind_data->variants.id_to_idx.find(vid);
				if (it == bind_data->variants.id_to_idx.end()) {
					throw InvalidInputException("read_pfile: variant '%s' not found in .pvar", vid);
				}
				bind_data->variant_indices.push_back(it->second);
			}
		} else {
			throw InvalidInputException("read_pfile: variants parameter must be LIST(VARCHAR) or LIST(INTEGER)");
		}

		bind_data->has_variant_filter = true;
	}

	// --- Build effective variant list (intersection of region + variant filter) ---
	{
		// Start with all variants if no filter, or variant_indices if variant filter active
		std::unordered_set<uint32_t> variant_set;
		bool use_variant_set = bind_data->has_variant_filter;
		if (use_variant_set) {
			for (auto idx : bind_data->variant_indices) {
				variant_set.insert(idx);
			}
		}

		bool any_filter_active = bind_data->region.active || bind_data->has_variant_filter;

		if (any_filter_active) {
			bind_data->has_effective_variant_list = true;

			for (uint32_t vidx = 0; vidx < bind_data->raw_variant_ct; vidx++) {
				// Check region filter
				if (bind_data->region.active) {
					if (bind_data->variants.chroms[vidx] != bind_data->region.chrom) {
						continue;
					}
					int64_t pos = bind_data->variants.positions[vidx];
					if (pos < bind_data->region.start || pos > bind_data->region.end) {
						continue;
					}
				}

				// Check variant filter
				if (use_variant_set && variant_set.find(vidx) == variant_set.end()) {
					continue;
				}

				bind_data->effective_variant_indices.push_back(vidx);
			}
		}
	}

	// --- Build output schema ---
	if (bind_data->tidy_mode) {
		// Tidy mode: variant columns + sample columns + genotype
		// Variant metadata
		names = {"CHROM", "POS", "ID", "REF", "ALT"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
		                LogicalType::VARCHAR};

		// Sample metadata columns from .psam
		auto &psam_header = bind_data->sample_metadata.header;
		for (idx_t i = 0; i < psam_header.column_names.size(); i++) {
			names.push_back(psam_header.column_names[i]);
			return_types.push_back(psam_header.column_types[i]);
			bind_data->tidy_sample_col_to_psam_col.push_back(i);
		}

		bind_data->tidy_sample_col_start = 5;

		// Genotype column (scalar, not list)
		names.push_back("genotype");
		return_types.push_back(LogicalType::TINYINT);
		bind_data->tidy_genotype_col = names.size() - 1;
		bind_data->tidy_total_cols = names.size();
	} else {
		// Default mode: same as read_pgen
		names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
		                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::TINYINT)};
	}

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PfileInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PfileBindData>();
	auto state = make_uniq<PfileGlobalState>();

	state->column_ids = input.column_ids;

	if (bind_data.tidy_mode) {
		state->total_effective_variants = bind_data.EffectiveVariantCt();
		// Check if genotype column is projected
		state->need_genotypes = false;
		for (auto col_id : input.column_ids) {
			if (col_id == bind_data.tidy_genotype_col) {
				state->need_genotypes = true;
				break;
			}
		}
	} else {
		state->total_variants = bind_data.EffectiveVariantCt();
		state->need_genotypes = false;
		for (auto col_id : input.column_ids) {
			if (col_id == PfileBindData::GENOTYPES_COL) {
				state->need_genotypes = true;
				break;
			}
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PfileInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                          GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PfileBindData>();
	auto &gstate = global_state->Cast<PfileGlobalState>();
	auto state = make_uniq<PfileLocalState>();

	if (!gstate.need_genotypes) {
		return std::move(state);
	}

	// --- Initialize per-thread PgenFileInfo + PgenReader ---
	plink2::PreinitPgfi(&state->pgfi);
	plink2::PreinitPgr(&state->pgr);

	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err =
	    plink2::PgfiInitPhase1(bind_data.pgen_path.c_str(), nullptr, bind_data.raw_variant_ct, bind_data.raw_sample_ct,
	                           &header_ctrl, &state->pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&state->pgfi, &cleanup_err);
		throw IOException("read_pfile: thread init failed (phase 1): %s", errstr_buf);
	}

	if (pgfi_alloc_cacheline_ct > 0) {
		state->pgfi_alloc_buf.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	}

	uint32_t max_vrec_width = 0;
	uintptr_t pgr_alloc_cacheline_ct = 0;

	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, state->pgfi.raw_variant_ct, &max_vrec_width, &state->pgfi,
	                             state->pgfi_alloc_buf.As<unsigned char>(), &pgr_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&state->pgfi, &cleanup_err);
		throw IOException("read_pfile: thread init failed (phase 2): %s", errstr_buf);
	}

	if (pgr_alloc_cacheline_ct > 0) {
		state->pgr_alloc_buf.Allocate(pgr_alloc_cacheline_ct * plink2::kCacheline);
	}

	err = plink2::PgrInit(bind_data.pgen_path.c_str(), max_vrec_width, &state->pgfi, &state->pgr,
	                      state->pgr_alloc_buf.As<unsigned char>());

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&state->pgr, &cleanup_err);
		cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&state->pgfi, &cleanup_err);
		throw IOException("read_pfile: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Allocate genovec buffer based on raw_sample_ct — pgenlib needs the full
	// genovec for internal decompression even when subsetting
	uint32_t genovec_sample_ct = bind_data.raw_sample_ct;
	uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(genovec_sample_ct);
	state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));
	std::memset(state->genovec_buf.ptr, 0, genovec_word_ct * sizeof(uintptr_t));

	state->genotype_bytes.resize(genovec_sample_ct);

	// Set up sample subsetting
	if (bind_data.has_sample_subset) {
		uintptr_t include_word_ct = plink2::DivUp(bind_data.raw_sample_ct, static_cast<uint32_t>(plink2::kBitsPerWord));
		state->sample_include_buf.Allocate(include_word_ct * sizeof(uintptr_t));
		auto *sample_include = state->sample_include_buf.As<uintptr_t>();
		std::memset(sample_include, 0, include_word_ct * sizeof(uintptr_t));

		for (auto idx : bind_data.sample_indices) {
			plink2::SetBit(idx, sample_include);
		}

		state->cumulative_popcounts_buf.Allocate(include_word_ct * sizeof(uint32_t));
		auto *cumulative_popcounts = state->cumulative_popcounts_buf.As<uint32_t>();
		plink2::FillCumulativePopcounts(sample_include, include_word_ct, cumulative_popcounts);

		plink2::PgrSetSampleSubsetIndex(cumulative_popcounts, &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Helper: resolve actual pgen variant index from effective position
// ---------------------------------------------------------------------------

static inline uint32_t ResolveVariantIdx(const PfileBindData &bind_data, uint32_t effective_pos) {
	if (bind_data.has_effective_variant_list) {
		return bind_data.effective_variant_indices[effective_pos];
	}
	return effective_pos;
}

// ---------------------------------------------------------------------------
// Scan: Default mode (variant-centric, same output as read_pgen)
// ---------------------------------------------------------------------------

static constexpr uint32_t PFILE_BATCH_SIZE = 128;

static void PfileDefaultScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();
	auto &gstate = data_p.global_state->Cast<PfileGlobalState>();
	auto &lstate = data_p.local_state->Cast<PfileLocalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t total_variants = gstate.total_variants;

	uint32_t output_sample_ct = bind_data.OutputSampleCt();

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(PFILE_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
		if (batch_start >= total_variants) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, total_variants);

		for (uint32_t effective_pos = batch_start; effective_pos < batch_end; effective_pos++) {
			uint32_t vidx = ResolveVariantIdx(bind_data, effective_pos);

			// Read genotype data if needed
			bool genotypes_read = false;
			if (gstate.need_genotypes && lstate.initialized) {
				const uintptr_t *sample_include =
				    bind_data.has_sample_subset ? lstate.sample_include_buf.As<uintptr_t>() : nullptr;

				plink2::PglErr err = plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
				                                    lstate.genovec_buf.As<uintptr_t>());

				if (err != plink2::kPglRetSuccess) {
					throw IOException("read_pfile: PgrGet failed for variant %u", vidx);
				}

				plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
				                             lstate.genotype_bytes.data());
				genotypes_read = true;
			}

			// Fill projected columns
			for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
				auto file_col = column_ids[out_col];
				if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
					continue;
				}

				auto &vec = output.data[out_col];

				switch (file_col) {
				case PfileBindData::CHROM_COL: {
					auto &val = bind_data.variants.chroms[vidx];
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case PfileBindData::POS_COL: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.positions[vidx];
					break;
				}
				case PfileBindData::ID_COL: {
					auto &val = bind_data.variants.ids[vidx];
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case PfileBindData::REF_COL: {
					auto &val = bind_data.variants.refs[vidx];
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case PfileBindData::ALT_COL: {
					auto &val = bind_data.variants.alts[vidx];
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case PfileBindData::GENOTYPES_COL: {
					if (!genotypes_read) {
						FlatVector::SetNull(vec, rows_emitted, true);
						break;
					}

					auto list_size = static_cast<idx_t>(output_sample_ct);
					auto current_offset = ListVector::GetListSize(vec);
					auto &entry = FlatVector::GetData<list_entry_t>(vec)[rows_emitted];
					entry.offset = current_offset;
					entry.length = list_size;

					ListVector::Reserve(vec, current_offset + list_size);
					auto &child = ListVector::GetEntry(vec);
					auto *child_data = FlatVector::GetData<int8_t>(child);
					auto &child_validity = FlatVector::Validity(child);

					for (idx_t s = 0; s < list_size; s++) {
						int8_t geno = lstate.genotype_bytes[s];
						if (geno == -9) {
							child_validity.SetInvalid(current_offset + s);
							child_data[current_offset + s] = 0;
						} else {
							child_data[current_offset + s] = geno;
						}
					}

					ListVector::SetListSize(vec, current_offset + list_size);
					break;
				}
				default:
					break;
				}
			}

			rows_emitted++;
		}
	}

	output.SetCardinality(rows_emitted);
}

// ---------------------------------------------------------------------------
// Scan: Tidy mode (one row per variant x sample)
// ---------------------------------------------------------------------------

//! Fill a sample metadata value into a DuckDB vector.
//! Handles SEX column (INTEGER), missing values, and parent columns.
static void FillSampleMetadataValue(Vector &vec, idx_t row_idx, const string &val, idx_t psam_col_idx,
                                    const PfileSampleMetadata &meta) {
	// SEX column: INTEGER type
	if (psam_col_idx == meta.sex_col_idx) {
		if (PfileIsMissingValue(val)) {
			FlatVector::SetNull(vec, row_idx, true);
		} else {
			try {
				int32_t sex_val = std::stoi(val);
				if (sex_val == 0) {
					FlatVector::SetNull(vec, row_idx, true);
				} else {
					FlatVector::GetData<int32_t>(vec)[row_idx] = sex_val;
				}
			} catch (...) {
				FlatVector::SetNull(vec, row_idx, true);
			}
		}
		return;
	}

	// PAT/MAT columns: "0" means unknown → NULL
	for (auto &parent_idx : meta.parent_col_indices) {
		if (psam_col_idx == parent_idx) {
			if (val == "0" || PfileIsMissingValue(val)) {
				FlatVector::SetNull(vec, row_idx, true);
			} else {
				FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, val);
			}
			return;
		}
	}

	// General VARCHAR columns
	if (PfileIsMissingValue(val)) {
		FlatVector::SetNull(vec, row_idx, true);
	} else {
		FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, val);
	}
}

static void PfileTidyScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();
	auto &gstate = data_p.global_state->Cast<PfileGlobalState>();
	auto &lstate = data_p.local_state->Cast<PfileLocalState>();

	if (lstate.tidy_done) {
		output.SetCardinality(0);
		return;
	}

	auto &column_ids = gstate.column_ids;
	uint32_t total_effective_variants = gstate.total_effective_variants;
	uint32_t output_sample_ct = bind_data.OutputSampleCt();

	// Build the list of sample indices to iterate in order
	// If subsetting, iterate over sample_indices in order.
	// If not subsetting, iterate 0..sample_ct-1.
	const vector<uint32_t> *sample_order = nullptr;
	vector<uint32_t> default_sample_order;
	if (bind_data.has_sample_subset) {
		sample_order = &bind_data.sample_indices;
	} else {
		default_sample_order.resize(output_sample_ct);
		for (uint32_t i = 0; i < output_sample_ct; i++) {
			default_sample_order[i] = i;
		}
		sample_order = &default_sample_order;
	}

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		if (lstate.tidy_current_variant_pos >= total_effective_variants) {
			lstate.tidy_done = true;
			break;
		}

		uint32_t vidx = ResolveVariantIdx(bind_data, lstate.tidy_current_variant_pos);

		// Load genotypes for current variant if not yet loaded
		if (!lstate.tidy_variant_loaded && gstate.need_genotypes && lstate.initialized) {
			const uintptr_t *sample_include =
			    bind_data.has_sample_subset ? lstate.sample_include_buf.As<uintptr_t>() : nullptr;

			plink2::PglErr err = plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
			                                    lstate.genovec_buf.As<uintptr_t>());

			if (err != plink2::kPglRetSuccess) {
				throw IOException("read_pfile: PgrGet failed for variant %u", vidx);
			}

			plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
			                             lstate.genotype_bytes.data());
			lstate.tidy_variant_loaded = true;
		}

		// Emit rows for samples within current variant
		while (lstate.tidy_current_sample < output_sample_ct && rows_emitted < STANDARD_VECTOR_SIZE) {

			uint32_t sample_file_idx = (*sample_order)[lstate.tidy_current_sample];

			// Fill projected columns
			for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
				auto file_col = column_ids[out_col];
				if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
					continue;
				}

				auto &vec = output.data[out_col];

				if (file_col < 5) {
					// Variant metadata columns
					switch (file_col) {
					case PfileBindData::CHROM_COL: {
						auto &val = bind_data.variants.chroms[vidx];
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						break;
					}
					case PfileBindData::POS_COL: {
						FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.positions[vidx];
						break;
					}
					case PfileBindData::ID_COL: {
						auto &val = bind_data.variants.ids[vidx];
						if (val.empty()) {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						}
						break;
					}
					case PfileBindData::REF_COL: {
						auto &val = bind_data.variants.refs[vidx];
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						break;
					}
					case PfileBindData::ALT_COL: {
						auto &val = bind_data.variants.alts[vidx];
						if (val.empty() || val == ".") {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						}
						break;
					}
					}
				} else if (file_col == bind_data.tidy_genotype_col) {
					// Genotype column (scalar TINYINT)
					if (gstate.need_genotypes && lstate.tidy_variant_loaded) {
						int8_t geno = lstate.genotype_bytes[lstate.tidy_current_sample];
						if (geno == -9) {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
						}
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
				} else if (file_col >= bind_data.tidy_sample_col_start && file_col < bind_data.tidy_genotype_col) {
					// Sample metadata column
					idx_t sample_col_rel = file_col - bind_data.tidy_sample_col_start;
					idx_t psam_col_idx = bind_data.tidy_sample_col_to_psam_col[sample_col_rel];

					if (sample_file_idx < bind_data.sample_metadata.rows.size()) {
						auto &sample_row = bind_data.sample_metadata.rows[sample_file_idx];
						if (psam_col_idx < sample_row.size()) {
							FillSampleMetadataValue(vec, rows_emitted, sample_row[psam_col_idx], psam_col_idx,
							                        bind_data.sample_metadata);
						} else {
							FlatVector::SetNull(vec, rows_emitted, true);
						}
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
				}
			}

			rows_emitted++;
			lstate.tidy_current_sample++;
		}

		// Advance to next variant if all samples emitted
		if (lstate.tidy_current_sample >= output_sample_ct) {
			lstate.tidy_current_variant_pos++;
			lstate.tidy_current_sample = 0;
			lstate.tidy_variant_loaded = false;
		}
	}

	output.SetCardinality(rows_emitted);
}

// ---------------------------------------------------------------------------
// Scan dispatch
// ---------------------------------------------------------------------------

static void PfileScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();

	if (bind_data.tidy_mode) {
		PfileTidyScan(context, data_p, output);
	} else {
		PfileDefaultScan(context, data_p, output);
	}
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPfileReader(ExtensionLoader &loader) {
	TableFunction read_pfile("read_pfile", {LogicalType::VARCHAR}, PfileScan, PfileBind, PfileInitGlobal,
	                         PfileInitLocal);

	read_pfile.projection_pushdown = true;

	read_pfile.named_parameters["pgen"] = LogicalType::VARCHAR;
	read_pfile.named_parameters["pvar"] = LogicalType::VARCHAR;
	read_pfile.named_parameters["psam"] = LogicalType::VARCHAR;
	read_pfile.named_parameters["tidy"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["dosages"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["phased"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["region"] = LogicalType::VARCHAR;
	// Accept ANY for samples and variants — type dispatch handled in bind
	read_pfile.named_parameters["samples"] = LogicalType::ANY;
	read_pfile.named_parameters["variants"] = LogicalType::ANY;

	loader.RegisterFunction(read_pfile);
}

} // namespace duckdb
