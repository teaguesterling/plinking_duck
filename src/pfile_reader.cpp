#include "pfile_reader.hpp"
#include "plink_common.hpp"
#include "plink_profile.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <pgenlib_read.h>
#include <pgenlib_ffi_support.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <unordered_set>

namespace duckdb {

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
// Sample metadata for genotype orient mode output
// ---------------------------------------------------------------------------

//! Extended sample info for genotype orient mode: includes all psam columns.
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

//! Load full sample metadata from a non-native source for genotype/sample orient modes.
//! Queries the source via Connection, populates both PfileSampleMetadata and SampleInfo.
static PfileSampleMetadata LoadPfileSampleMetadataFromSource(ClientContext &context, const string &source,
                                                             SampleInfo &sample_info_out) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto escaped = StringUtil::Replace(source, "'", "''");
	auto result = conn.Query("SELECT * FROM '" + escaped + "'");
	if (result->HasError()) {
		throw IOException("read_pfile: failed to read source '%s': %s", source, result->GetError());
	}

	PfileSampleMetadata meta;

	// Build header from result schema
	meta.header.format = PsamFormat::PSAM_IID;
	meta.header.column_names = result->names;
	meta.header.column_types = result->types;

	// Find special column indices (case-insensitive)
	idx_t iid_idx = DConstants::INVALID_INDEX;
	idx_t fid_idx = DConstants::INVALID_INDEX;
	for (idx_t i = 0; i < result->names.size(); i++) {
		auto lower = StringUtil::Lower(result->names[i]);
		if (lower == "sex") {
			meta.sex_col_idx = i;
		} else if (lower == "pat" || lower == "mat") {
			meta.parent_col_indices.push_back(i);
		}
		if (lower == "iid") {
			iid_idx = i;
		} else if (lower == "fid") {
			fid_idx = i;
		}
	}

	if (iid_idx == DConstants::INVALID_INDEX) {
		throw IOException("read_pfile: source '%s' has no IID column (found: %s)", source,
		                  StringUtil::Join(result->names, ", "));
	}

	bool has_fid = (fid_idx != DConstants::INVALID_INDEX);

	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		for (idx_t row = 0; row < chunk->size(); row++) {
			vector<string> fields;
			for (idx_t col = 0; col < chunk->ColumnCount(); col++) {
				auto val = chunk->GetValue(col, row);
				fields.push_back(val.IsNull() ? "" : val.ToString());
			}

			// iid_to_idx is lazy (see SampleInfo::EnsureIidMap).
			sample_info_out.iids.push_back(fields[iid_idx]);
			if (has_fid) {
				sample_info_out.fids.push_back(fields[fid_idx]);
			}

			meta.rows.push_back(std::move(fields));
		}
	}

	sample_info_out.sample_ct = sample_info_out.iids.size();
	return meta;
}

//! Load full sample metadata including all columns for genotype orient mode output.
//! Also populates sample_info to avoid a separate file read for LoadSampleInfo.
static PfileSampleMetadata LoadPfileSampleMetadata(ClientContext &context, const string &path,
                                                   SampleInfo &sample_info_out) {
	// Non-native sources: dispatch to source loader
	if (!IsNativePlinkFormat(path)) {
		return LoadPfileSampleMetadataFromSource(context, path, sample_info_out);
	}

	auto lines = ReadFileLines(context, path);

	PfileSampleMetadata meta;

	if (lines.empty()) {
		throw IOException("read_pfile: .psam/.fam file '%s' is empty", path);
	}

	// Parse header from first line (avoids re-reading the file)
	meta.header = ParsePsamHeader(context, path);

	// Find special column indices
	idx_t iid_idx = DConstants::INVALID_INDEX;
	idx_t fid_idx = DConstants::INVALID_INDEX;
	for (idx_t i = 0; i < meta.header.column_names.size(); i++) {
		auto &name = meta.header.column_names[i];
		if (name == "SEX") {
			meta.sex_col_idx = i;
		} else if (name == "PAT" || name == "MAT") {
			meta.parent_col_indices.push_back(i);
		}
		if (name == "IID") {
			iid_idx = i;
		} else if (name == "FID") {
			fid_idx = i;
		}
	}

	if (iid_idx == DConstants::INVALID_INDEX) {
		throw IOException("read_pfile: .psam/.fam file '%s' has no IID column", path);
	}

	bool has_fid = (fid_idx != DConstants::INVALID_INDEX);

	// Skip header/comment lines: .fam has no header (start at 0),
	// .psam may have ## comment lines before the # header line.
	idx_t data_start = 0;
	if (meta.header.format != PsamFormat::FAM) {
		while (data_start < lines.size() && !lines[data_start].empty() && lines[data_start][0] == '#') {
			data_start++;
		}
	}

	for (idx_t i = data_start; i < lines.size(); i++) {
		auto &line = lines[i];
		if (line.empty()) {
			continue;
		}

		vector<string> fields;
		if (meta.header.format == PsamFormat::FAM) {
			fields = SplitWhitespaceLine(line);
		} else {
			fields = SplitTabLine(line);
		}

		// iid_to_idx is built lazily on demand (see SampleInfo::EnsureIidMap).
		if (iid_idx < fields.size()) {
			sample_info_out.iids.push_back(fields[iid_idx]);
			if (has_fid && fid_idx < fields.size()) {
				sample_info_out.fids.push_back(fields[fid_idx]);
			}
		}

		meta.rows.push_back(std::move(fields));
	}

	sample_info_out.sample_ct = sample_info_out.iids.size();
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
	VariantMetadataIndex variants;

	// Sample metadata (basic: for IID lookups and default mode)
	SampleInfo sample_info;

	// Full sample metadata (for genotype orient mode column output)
	PfileSampleMetadata sample_metadata;

	// pgenlib header info
	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Mode
	OrientMode orient_mode = OrientMode::VARIANT;

	// Options
	bool include_dosages = false;
	bool include_phased = false;
	GenotypeMode genotype_mode = GenotypeMode::ARRAY;

	// Sample subsetting
	bool has_sample_subset = false;
	vector<uint32_t> sample_indices; // 0-based indices into .pgen sample order
	uint32_t subset_sample_ct = 0;

	// Region filtering
	RegionFilter region;

	// Variant filtering
	bool has_variant_filter = false;
	vector<uint32_t> variant_indices; // sorted indices of variants to include

	// Count-based filtering (af_range, ac_range)
	CountFilter count_filter;
	unique_ptr<SampleSubset> count_filter_subset; // shared bind-time subset for PgrGetCounts

	// Genotype range filtering (genotype_range)
	GenotypeRangeFilter genotype_filter;

	// Effective variant list (after region + variant filter intersection)
	// If empty and no filters active, scan all variants sequentially.
	// If non-empty, these are the specific variant indices to scan.
	bool has_effective_variant_list = false;
	vector<uint32_t> effective_variant_indices;

	// Sample-orient mode: pre-read genotype matrix (variant × sample)
	// genotype_matrix[effective_vidx][sample_idx] = genotype value (-9 = missing)
	vector<vector<int8_t>> genotype_matrix;
	// Dosage variant of genotype_matrix (-9.0 = missing)
	vector<vector<double>> dosage_matrix;

	// Sample-orient column layout
	// Sample metadata columns start at index 0, genotypes column is last
	idx_t sample_orient_genotypes_col = 0; // computed in bind (ARRAY/LIST only)
	idx_t sample_orient_total_cols = 0;    // computed in bind
	// Mapping from sample-orient output column index to psam file column index
	vector<idx_t> sample_orient_col_to_psam_col;

	// Columns mode layout (genotypes := 'columns')
	vector<string> genotype_column_names;     // IIDs (variant orient) or variant IDs (sample orient)
	idx_t columns_mode_first_geno_col = 0;    // first genotype column index in schema
	uint32_t columns_mode_geno_col_count = 0; // number of genotype columns

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
	idx_t geno_sample_col_start = 5;
	idx_t geno_genotype_col = 0; // computed in bind
	idx_t geno_total_cols = 0;   // computed in bind

	// Mapping from genotype orient output sample column index (relative to geno_sample_col_start)
	// to psam file column index
	vector<idx_t> geno_sample_col_to_psam_col;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PfileGlobalState : public GlobalTableFunctionState {
	// For variant/sample-orient mode: atomic counter
	std::atomic<uint32_t> next_idx {0};
	uint32_t total_count = 0; // variants (variant/genotype mode) or samples (sample mode)

	// Projection
	bool need_genotypes = false;
	bool need_pgen_reader = false; // genotypes OR count filter
	vector<column_t> column_ids;
	OrientMode orient_mode = OrientMode::VARIANT;
	uint32_t max_threads_config = 0;

	idx_t MaxThreads() const override {
		if (orient_mode == OrientMode::GENOTYPE) {
			// Each variant fans out to N sample rows — use smaller batch size (64)
			// for better load balancing when filters skip variants unevenly
			return ApplyMaxThreadsCap(total_count / 64 + 1, max_threads_config);
		}
		// Variant and sample modes support parallel scan
		return ApplyMaxThreadsCap(total_count / 1000 + 1, max_threads_config);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PfileLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;
	AlignedBuffer genovec_buf;
	AlignedBuffer sample_include_buf;
	AlignedBuffer cumulative_popcounts_buf;

	vector<int8_t> genotype_bytes;

	// Phase buffers (for phased := true)
	AlignedBuffer phasepresent_buf;
	AlignedBuffer phaseinfo_buf;
	vector<int8_t> phased_pairs;
	uint32_t phasepresent_ct = 0;

	// Dosage buffers (for dosages := true)
	AlignedBuffer dosage_present_buf;
	AlignedBuffer dosage_main_buf;
	vector<double> dosage_doubles;

	plink2::PgrSampleSubsetIndex pssi;

	bool initialized = false;

	// Genotype-orient batch claiming state
	uint32_t batch_start = 0;               // first effective variant index in current batch
	uint32_t batch_end = 0;                 // one-past-end effective variant index in current batch
	uint32_t current_variant_in_batch = 0;  // current effective variant index within batch
	uint32_t current_sample_in_variant = 0; // sample index within current variant
	bool batch_variant_loaded = false;      // genotypes decoded for current variant in batch
	bool geno_range_all_pass = true;        // per-variant flag for genotype_range optimization
	bool batch_exhausted = true;            // true initially to trigger first batch claim

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
	BindPhaseTimer bind_timer("PfileBind(total)");
	auto bind_data = make_uniq<PfileBindData>();

	// --- Resolve file paths ---
	// First positional argument is the prefix (or empty if only named params)
	string prefix;
	if (!input.inputs.empty()) {
		prefix = input.inputs[0].GetValue<string>();
	}

	auto &fs = FileSystem::GetFileSystem(context);

	// Named parameters can override individual paths
	string orient_str;
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pgen") {
			bind_data->pgen_path = kv.second.GetValue<string>();
		} else if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "orient") {
			orient_str = kv.second.GetValue<string>();
		} else if (kv.first == "dosages") {
			bind_data->include_dosages = kv.second.GetValue<bool>();
		} else if (kv.first == "phased") {
			bind_data->include_phased = kv.second.GetValue<bool>();
		} else if (kv.first == "region") {
			bind_data->region = ParseRegion(kv.second.GetValue<string>());
		}
		// samples, variants, genotypes, af_range, ac_range handled after pgenlib init
	}

	bind_data->orient_mode = ResolveOrientMode(orient_str, "read_pfile");

	if (bind_data->include_dosages && bind_data->include_phased) {
		throw InvalidInputException("read_pfile: dosages and phased cannot both be true");
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
			bind_data->pvar_path = FindCompanionFileWithParquet(context, fs, prefix + ".pgen", {".pvar", ".bim"});
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
			bind_data->pvar_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".pvar", ".bim"});
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
			bind_data->psam_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".psam", ".fam"});
		}
		if (bind_data->psam_path.empty()) {
			throw InvalidInputException("read_pfile: cannot find .psam or .fam file for '%s' "
			                            "(use psam := 'path' to specify explicitly)",
			                            prefix.empty() ? bind_data->pgen_path : prefix);
		}
	}

	// --- Initialize pgenlib (Phase 1) to get counts ---
	bind_timer.Note("file paths resolved, starting pgenlib init");
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
	AlignedBuffer pgfi_alloc;
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
	bind_timer.Note("pgenlib init done (raw_variant_ct=%u, raw_sample_ct=%u), loading variant metadata",
	                bind_data->raw_variant_ct, bind_data->raw_sample_ct);
	// Parquet + region: push the WHERE into the parquet scan so we only materialize
	// the region's variants instead of all N (huge win at 170M).
	if (bind_data->region.active && IsParquetFile(bind_data->pvar_path)) {
		idx_t total_ct = GetParquetRowCount(context, bind_data->pvar_path);
		bind_data->variants = LoadVariantMetadataFromParquetRegion(context, bind_data->pvar_path,
		                                                           bind_data->region.chrom, bind_data->region.start,
		                                                           bind_data->region.end, total_ct, "read_pfile");
	} else {
		bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "read_pfile");
	}
	bind_timer.Note("variant metadata loaded (%llu total variants, %llu loaded)",
	                (unsigned long long)bind_data->variants.variant_ct,
	                (unsigned long long)bind_data->variants.chroms.size());

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("read_pfile: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info ---
	// Determine whether the query actually needs IID strings. If not, we only
	// need the row count — at biobank scale this saves ~600ms of string copies.
	// Need IIDs when:
	//   * orient = genotype/sample (full per-sample metadata)
	//   * samples filter contains VARCHAR (IID-string lookup)
	//   * genotypes := 'columns' or 'struct' in variant orient (column names from IIDs)
	bool needs_iids = false;
	if (bind_data->orient_mode == OrientMode::GENOTYPE || bind_data->orient_mode == OrientMode::SAMPLE) {
		needs_iids = true;
	}
	{
		auto it = input.named_parameters.find("samples");
		if (it != input.named_parameters.end()) {
			auto &child_type = ListType::GetChildType(it->second.type());
			if (child_type.id() == LogicalTypeId::VARCHAR) {
				needs_iids = true;
			}
		}
	}
	{
		auto it = input.named_parameters.find("genotypes");
		if (it != input.named_parameters.end() && bind_data->orient_mode == OrientMode::VARIANT) {
			auto gv = StringUtil::Lower(it->second.GetValue<string>());
			if (gv == "columns" || gv == "struct") {
				needs_iids = true;
			}
		}
	}

	bind_timer.Note("loading sample metadata (needs_iids=%s)", needs_iids ? "yes" : "no");
	if (bind_data->orient_mode == OrientMode::GENOTYPE || bind_data->orient_mode == OrientMode::SAMPLE) {
		bind_data->sample_metadata = LoadPfileSampleMetadata(context, bind_data->psam_path, bind_data->sample_info);
	} else if (needs_iids) {
		bind_data->sample_info = LoadSampleMetadata(context, bind_data->psam_path);
	} else {
		bind_data->sample_info = LoadSampleCount(context, bind_data->psam_path);
	}
	bind_timer.Note("sample metadata loaded (%llu samples, iids=%llu)",
	                (unsigned long long)bind_data->sample_info.sample_ct,
	                (unsigned long long)bind_data->sample_info.iids.size());

	if (bind_data->sample_info.sample_ct != bind_data->raw_sample_ct) {
		throw InvalidInputException("read_pfile: sample count mismatch: .pgen has %u samples, "
		                            ".psam/.fam '%s' has %llu samples",
		                            bind_data->raw_sample_ct, bind_data->psam_path,
		                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
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
			// Lazily build iid_to_idx only when VARCHAR filter is actually used.
			bind_data->sample_info.EnsureIidMap("read_pfile: " + bind_data->psam_path);
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
		// genotype_bytes[i] corresponds to sample_indices[i] in genotype orient mode.
		std::sort(bind_data->sample_indices.begin(), bind_data->sample_indices.end());

		bind_data->has_sample_subset = true;
		bind_data->subset_sample_ct = static_cast<uint32_t>(bind_data->sample_indices.size());
	}

	// --- Process variants parameter ---
	auto variants_it = input.named_parameters.find("variants");
	if (variants_it != input.named_parameters.end()) {
		bind_data->variant_indices =
		    ResolveVariantsParameter(variants_it->second, bind_data->variants, bind_data->raw_variant_ct, "read_pfile");
		bind_data->has_variant_filter = true;
	}

	// --- Build effective variant list (intersection of region + variant filter) ---
	bind_timer.Note("building effective variant list");
	{
		std::unordered_set<uint32_t> variant_set;
		bool use_variant_set = bind_data->has_variant_filter;
		if (use_variant_set) {
			for (auto idx : bind_data->variant_indices) {
				variant_set.insert(idx);
			}
		}

		bool any_filter_active = bind_data->region.active || bind_data->has_variant_filter;

		if (any_filter_active) {
			BindPhaseTimer evl_timer("effective-variant-list-scan");
			bind_data->has_effective_variant_list = true;

			// Case A: sparse pvar index (parquet region pushdown). The loaded
			// subset IS the region-matched set; iterate vidx_map keys directly.
			if (!bind_data->variants.vidx_map.empty()) {
				bind_data->effective_variant_indices.reserve(bind_data->variants.vidx_map.size());
				for (auto &kv : bind_data->variants.vidx_map) {
					if (use_variant_set && variant_set.find(kv.first) == variant_set.end()) {
						continue;
					}
					bind_data->effective_variant_indices.push_back(kv.first);
				}
				std::sort(bind_data->effective_variant_indices.begin(), bind_data->effective_variant_indices.end());
				evl_timer.Note("sparse pvar: %llu region variants pre-filtered; %llu passed",
				               (unsigned long long)bind_data->variants.vidx_map.size(),
				               (unsigned long long)bind_data->effective_variant_indices.size());
			} else if (bind_data->region.active && !bind_data->variants.chrom_offsets.empty()) {
				// Case B: dense + region + chrom_offsets → O(log N) binary-search bounds.
				auto it = bind_data->variants.chrom_offsets.find(bind_data->region.chrom);
				if (it != bind_data->variants.chrom_offsets.end()) {
					idx_t lo_local = it->second.first;
					idx_t hi_local = it->second.second;
					auto &positions = bind_data->variants.positions;
					// binary search for pos >= region.start
					idx_t lo = lo_local, hi = hi_local;
					while (lo < hi) {
						idx_t mid = lo + (hi - lo) / 2;
						if (positions[mid] < bind_data->region.start) {
							lo = mid + 1;
						} else {
							hi = mid;
						}
					}
					idx_t start_idx = lo;
					// binary search for pos > region.end
					lo = lo_local;
					hi = hi_local;
					while (lo < hi) {
						idx_t mid = lo + (hi - lo) / 2;
						if (positions[mid] <= bind_data->region.end) {
							lo = mid + 1;
						} else {
							hi = mid;
						}
					}
					idx_t end_idx = lo;
					bind_data->effective_variant_indices.reserve(end_idx - start_idx);
					for (idx_t vidx = start_idx; vidx < end_idx; vidx++) {
						if (use_variant_set && variant_set.find(static_cast<uint32_t>(vidx)) == variant_set.end()) {
							continue;
						}
						bind_data->effective_variant_indices.push_back(static_cast<uint32_t>(vidx));
					}
				}
				evl_timer.Note("dense+region: chrom_offsets+bsearch → %llu passed",
				               (unsigned long long)bind_data->effective_variant_indices.size());
			} else {
				// Case C: variant filter only (no region), or region without chrom_offsets.
				// Full scan with cheap columnar access.
				for (uint32_t vidx = 0; vidx < bind_data->raw_variant_ct; vidx++) {
					if (bind_data->region.active) {
						if (bind_data->variants.GetChrom(vidx) != bind_data->region.chrom) {
							continue;
						}
						int64_t pos = bind_data->variants.GetPos(vidx);
						if (pos < bind_data->region.start || pos > bind_data->region.end) {
							continue;
						}
					}
					if (use_variant_set && variant_set.find(vidx) == variant_set.end()) {
						continue;
					}
					bind_data->effective_variant_indices.push_back(vidx);
				}
				evl_timer.Note("linear fallback: %u variants scanned, %llu passed", bind_data->raw_variant_ct,
				               (unsigned long long)bind_data->effective_variant_indices.size());
			}
		}
	}

	// --- Parse count filters (af_range, ac_range) ---
	{
		auto af_it = input.named_parameters.find("af_range");
		if (af_it != input.named_parameters.end()) {
			bind_data->count_filter.af_filter = ParseRangeFilter(af_it->second, "af_range", 0.0, 1.0, "read_pfile");
		}
		auto ac_it = input.named_parameters.find("ac_range");
		if (ac_it != input.named_parameters.end()) {
			bind_data->count_filter.ac_filter = ParseRangeFilter(
			    ac_it->second, "ac_range", 0.0, static_cast<double>(2 * bind_data->OutputSampleCt()), "read_pfile");
		}
	}

	// --- Parse genotype_range filter ---
	{
		auto gr_it = input.named_parameters.find("genotype_range");
		if (gr_it != input.named_parameters.end()) {
			if (bind_data->include_dosages) {
				throw InvalidInputException("read_pfile: genotype_range is incompatible with dosages := true");
			}
			bind_data->genotype_filter.range =
			    ParseRangeFilter(gr_it->second, "genotype_range", 0.0, 2.0, "read_pfile");
			bind_data->genotype_filter.active = bind_data->genotype_filter.range.active;
		}
	}

	// Build shared SampleSubset for PgrGetCounts if needed.
	// Required when af_range/ac_range or genotype_range is active with sample subsetting,
	// because PgrGetCounts needs sample_include + interleaved_vec for correct counting.
	// Must come after both count_filter and genotype_filter are parsed.
	if ((bind_data->count_filter.HasFilter() || bind_data->genotype_filter.active) && bind_data->has_sample_subset &&
	    !bind_data->count_filter_subset) {
		bind_data->count_filter_subset =
		    make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, bind_data->sample_indices));
	}

	// --- Build output schema ---
	if (bind_data->orient_mode == OrientMode::GENOTYPE) {
		// Check for incompatible genotypes modes before building schema
		auto genotypes_it_geno = input.named_parameters.find("genotypes");
		if (genotypes_it_geno != input.named_parameters.end()) {
			auto gval = StringUtil::Lower(genotypes_it_geno->second.GetValue<string>());
			if (gval == "columns") {
				throw InvalidInputException(
				    "read_pfile: genotypes := 'columns' is not compatible with orient := 'genotype' "
				    "(genotype mode already produces scalar output)");
			}
			if (gval == "struct") {
				throw InvalidInputException(
				    "read_pfile: genotypes := 'struct' is not compatible with orient := 'genotype' "
				    "(genotype mode already produces scalar output)");
			}
			if (gval == "counts" || gval == "stats") {
				throw InvalidInputException("read_pfile: genotypes := '%s' is not compatible with orient := 'genotype' "
				                            "(aggregate modes require orient := 'variant' or 'sample')",
				                            gval);
			}
		}

		// Genotype mode: variant columns + sample columns + genotype
		// Variant metadata
		names = {"CHROM", "POS", "ID", "REF", "ALT"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
		                LogicalType::VARCHAR};

		// Sample metadata columns from .psam
		auto &psam_header = bind_data->sample_metadata.header;
		for (idx_t i = 0; i < psam_header.column_names.size(); i++) {
			names.push_back(psam_header.column_names[i]);
			return_types.push_back(psam_header.column_types[i]);
			bind_data->geno_sample_col_to_psam_col.push_back(i);
		}

		bind_data->geno_sample_col_start = 5;

		// Genotype column (scalar TINYINT, ARRAY(TINYINT,2) when phased, or DOUBLE when dosages)
		names.push_back("genotype");
		return_types.push_back(bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
		                       : bind_data->include_dosages ? LogicalType::DOUBLE
		                                                    : LogicalType::TINYINT);
		bind_data->geno_genotype_col = names.size() - 1;
		bind_data->geno_total_cols = names.size();
	} else if (bind_data->orient_mode == OrientMode::SAMPLE) {
		// Apply count filter in sample-orient mode before building schema,
		// so ARRAY dimension reflects filtered variant count.
		vector<bool> genotype_range_all_pass;
		if (bind_data->count_filter.HasFilter() || bind_data->genotype_filter.active) {
			// Init temporary PgenReader for count filtering
			plink2::PgenFileInfo cf_pgfi;
			plink2::PreinitPgfi(&cf_pgfi);

			char cf_errstr[plink2::kPglErrstrBufBlen];
			plink2::PgenHeaderCtrl cf_hdr;
			uintptr_t cf_pgfi_alloc_ct = 0;

			err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr, bind_data->raw_variant_ct,
			                             bind_data->raw_sample_ct, &cf_hdr, &cf_pgfi, &cf_pgfi_alloc_ct, cf_errstr);
			if (err != plink2::kPglRetSuccess) {
				plink2::PglErr ce = plink2::kPglRetSuccess;
				plink2::CleanupPgfi(&cf_pgfi, &ce);
				throw IOException("read_pfile: failed to init PgenReader for count filter: %s", cf_errstr);
			}

			AlignedBuffer cf_pgfi_alloc;
			if (cf_pgfi_alloc_ct > 0) {
				cf_pgfi_alloc.Allocate(cf_pgfi_alloc_ct * plink2::kCacheline);
			}

			uint32_t cf_max_vrec = 0;
			uintptr_t cf_pgr_alloc_ct = 0;
			err = plink2::PgfiInitPhase2(cf_hdr, 0, 0, 0, 0, cf_pgfi.raw_variant_ct, &cf_max_vrec, &cf_pgfi,
			                             cf_pgfi_alloc.As<unsigned char>(), &cf_pgr_alloc_ct, cf_errstr);
			if (err != plink2::kPglRetSuccess) {
				plink2::PglErr ce = plink2::kPglRetSuccess;
				plink2::CleanupPgfi(&cf_pgfi, &ce);
				throw IOException("read_pfile: count filter PgenReader init failed (phase 2): %s", cf_errstr);
			}

			plink2::PgenReader cf_pgr;
			plink2::PreinitPgr(&cf_pgr);
			AlignedBuffer cf_pgr_alloc;
			if (cf_pgr_alloc_ct > 0) {
				cf_pgr_alloc.Allocate(cf_pgr_alloc_ct * plink2::kCacheline);
			}

			err = plink2::PgrInit(bind_data->pgen_path.c_str(), cf_max_vrec, &cf_pgfi, &cf_pgr,
			                      cf_pgr_alloc.As<unsigned char>());
			if (err != plink2::kPglRetSuccess) {
				plink2::PglErr ce = plink2::kPglRetSuccess;
				plink2::CleanupPgr(&cf_pgr, &ce);
				ce = plink2::kPglRetSuccess;
				plink2::CleanupPgfi(&cf_pgfi, &ce);
				throw IOException("read_pfile: PgrInit failed for count filter");
			}

			// Set up sample subsetting for PgrGetCounts
			plink2::PgrSampleSubsetIndex cf_pssi;
			if (bind_data->has_sample_subset && bind_data->count_filter_subset) {
				plink2::PgrSetSampleSubsetIndex(bind_data->count_filter_subset->CumulativePopcounts(), &cf_pgr,
				                                &cf_pssi);
			} else {
				plink2::PgrClearSampleSubsetIndex(&cf_pgr, &cf_pssi);
			}

			const uintptr_t *cf_si = (bind_data->has_sample_subset && bind_data->count_filter_subset)
			                             ? bind_data->count_filter_subset->SampleInclude()
			                             : nullptr;
			const uintptr_t *cf_iv = (bind_data->has_sample_subset && bind_data->count_filter_subset)
			                             ? bind_data->count_filter_subset->InterleavedVec()
			                             : nullptr;
			uint32_t cf_sc = bind_data->has_sample_subset ? bind_data->subset_sample_ct : bind_data->raw_sample_ct;

			// Iterate effective variants and filter by count
			uint32_t pre_filter_ct = bind_data->EffectiveVariantCt();
			vector<uint32_t> filtered_indices;
			filtered_indices.reserve(pre_filter_ct);

			for (uint32_t ev = 0; ev < pre_filter_ct; ev++) {
				uint32_t vidx = bind_data->has_effective_variant_list ? bind_data->effective_variant_indices[ev] : ev;

				STD_ARRAY_DECL(uint32_t, 4, genocounts);
				plink2::PglErr cf_err = plink2::PgrGetCounts(cf_si, cf_iv, cf_pssi, cf_sc, vidx, &cf_pgr, genocounts);
				if (cf_err != plink2::kPglRetSuccess) {
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&cf_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&cf_pgfi, &ce);
					throw IOException("read_pfile: PgrGetCounts failed for variant %u during count filter", vidx);
				}

				auto pf = CheckPreDecompFilters(bind_data->count_filter, bind_data->genotype_filter, genocounts, cf_sc);
				if (pf.skip) {
					continue;
				}
				bool all_pass = pf.all_pass;

				filtered_indices.push_back(vidx);
				genotype_range_all_pass.push_back(all_pass);
			}

			bind_data->effective_variant_indices = std::move(filtered_indices);
			bind_data->has_effective_variant_list = true;

			// Cleanup
			{
				plink2::PglErr ce = plink2::kPglRetSuccess;
				plink2::CleanupPgr(&cf_pgr, &ce);
				ce = plink2::kPglRetSuccess;
				plink2::CleanupPgfi(&cf_pgfi, &ce);
			}
		}

		// Sample-orient mode: sample metadata columns + genotypes array/list
		auto &psam_header = bind_data->sample_metadata.header;
		for (idx_t i = 0; i < psam_header.column_names.size(); i++) {
			names.push_back(psam_header.column_names[i]);
			return_types.push_back(psam_header.column_types[i]);
			bind_data->sample_orient_col_to_psam_col.push_back(i);
		}

		// Resolve genotype output mode — dimension is effective_variant_ct
		uint32_t effective_variant_ct = bind_data->EffectiveVariantCt();
		string genotypes_str = "auto";
		auto genotypes_it = input.named_parameters.find("genotypes");
		if (genotypes_it != input.named_parameters.end()) {
			genotypes_str = genotypes_it->second.GetValue<string>();
		}
		bind_data->genotype_mode = ResolveGenotypeMode(genotypes_str, effective_variant_ct, "read_pfile");

		// Validate incompatible combinations for aggregate modes
		if (IsAggregateGenotypeMode(bind_data->genotype_mode)) {
			if (bind_data->include_phased) {
				throw InvalidInputException("read_pfile: genotypes := '%s' is incompatible with phased := true",
				                            bind_data->genotype_mode == GenotypeMode::COUNTS ? "counts" : "stats");
			}
			if (bind_data->include_dosages) {
				throw InvalidInputException("read_pfile: genotypes := '%s' is incompatible with dosages := true",
				                            bind_data->genotype_mode == GenotypeMode::COUNTS ? "counts" : "stats");
			}
		}

		if (bind_data->genotype_mode == GenotypeMode::COLUMNS) {
			// Columns mode: one scalar TINYINT column per effective variant
			bind_data->columns_mode_first_geno_col = names.size();
			bind_data->columns_mode_geno_col_count = effective_variant_ct;

			// Check for duplicate variant IDs (variant IDs can duplicate unlike sample IIDs)
			std::unordered_set<string> seen_names;
			for (uint32_t ev = 0; ev < effective_variant_ct; ev++) {
				uint32_t vidx = bind_data->has_effective_variant_list ? bind_data->effective_variant_indices[ev] : ev;
				auto id = bind_data->variants.GetId(vidx);
				string col_name;
				if (id.empty()) {
					col_name =
					    bind_data->variants.GetChrom(vidx) + ":" + std::to_string(bind_data->variants.GetPos(vidx));
				} else {
					col_name = id;
				}
				if (!seen_names.insert(col_name).second) {
					throw InvalidInputException(
					    "read_pfile: genotypes := 'columns' with orient := 'sample' requires unique variant "
					    "identifiers, but '%s' appears more than once. Use variants := [...] to select unique "
					    "variants.",
					    col_name);
				}
				bind_data->genotype_column_names.push_back(col_name);
				names.push_back(col_name);
				return_types.push_back(bind_data->include_dosages ? LogicalType::DOUBLE : LogicalType::TINYINT);
			}

			bind_data->sample_orient_total_cols = names.size();
		} else if (bind_data->genotype_mode == GenotypeMode::STRUCT) {
			// STRUCT mode: single genotypes column with one field per effective variant
			child_list_t<LogicalType> struct_children;
			LogicalType field_type = bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
			                         : bind_data->include_dosages ? LogicalType::DOUBLE
			                                                      : LogicalType::TINYINT;

			std::unordered_set<string> seen_names;
			for (uint32_t ev = 0; ev < effective_variant_ct; ev++) {
				uint32_t vidx = bind_data->has_effective_variant_list ? bind_data->effective_variant_indices[ev] : ev;
				auto id = bind_data->variants.GetId(vidx);
				string col_name;
				if (id.empty()) {
					col_name =
					    bind_data->variants.GetChrom(vidx) + ":" + std::to_string(bind_data->variants.GetPos(vidx));
				} else {
					col_name = id;
				}
				if (!seen_names.insert(col_name).second) {
					throw InvalidInputException(
					    "read_pfile: genotypes := 'struct' with orient := 'sample' requires unique variant "
					    "identifiers, but '%s' appears more than once. Use variants := [...] to select unique "
					    "variants.",
					    col_name);
				}
				bind_data->genotype_column_names.push_back(col_name);
				struct_children.push_back({col_name, field_type});
			}

			names.push_back("genotypes");
			return_types.push_back(LogicalType::STRUCT(std::move(struct_children)));
			bind_data->sample_orient_genotypes_col = names.size() - 1;
			bind_data->sample_orient_total_cols = names.size();
		} else if (bind_data->genotype_mode == GenotypeMode::COUNTS) {
			names.push_back("genotypes");
			return_types.push_back(MakeGenotypeCountsType());
			bind_data->sample_orient_genotypes_col = names.size() - 1;
			bind_data->sample_orient_total_cols = names.size();
		} else if (bind_data->genotype_mode == GenotypeMode::STATS) {
			names.push_back("genotypes");
			return_types.push_back(MakeGenotypeStatsType());
			bind_data->sample_orient_genotypes_col = names.size() - 1;
			bind_data->sample_orient_total_cols = names.size();
		} else {
			LogicalType sample_elem_type = bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
			                               : bind_data->include_dosages ? LogicalType::DOUBLE
			                                                            : LogicalType::TINYINT;
			LogicalType geno_type = bind_data->genotype_mode == GenotypeMode::ARRAY
			                            ? LogicalType::ARRAY(sample_elem_type, effective_variant_ct)
			                            : LogicalType::LIST(sample_elem_type);

			names.push_back("genotypes");
			return_types.push_back(geno_type);

			bind_data->sample_orient_genotypes_col = names.size() - 1;
			bind_data->sample_orient_total_cols = names.size();
		}

		// --- Pre-read all genotypes into matrix ---
		uint32_t output_sample_ct = bind_data->OutputSampleCt();
		uint64_t matrix_size = static_cast<uint64_t>(effective_variant_ct) * static_cast<uint64_t>(output_sample_ct);

		// Check configurable matrix size limit
		int64_t max_elements = 16LL * 1024 * 1024 * 1024; // default 16G
		Value max_elements_val;
		if (context.TryGetCurrentSetting("plinking_max_matrix_elements", max_elements_val)) {
			max_elements = max_elements_val.GetValue<int64_t>();
		}
		if (matrix_size > static_cast<uint64_t>(max_elements)) {
			throw InvalidInputException("read_pfile: orient := 'sample' would require %llu genotype values "
			                            "(%u variants x %u samples, limit: %lld). "
			                            "Use variants := [...] or samples := [...] to reduce, "
			                            "or SET plinking_max_matrix_elements = <higher value>.",
			                            static_cast<unsigned long long>(matrix_size), effective_variant_ct,
			                            output_sample_ct, static_cast<long long>(max_elements));
		}

		// Init temporary PgenReader for pre-reading
		plink2::PgenFileInfo tmp_pgfi;
		plink2::PreinitPgfi(&tmp_pgfi);

		char errstr_buf2[plink2::kPglErrstrBufBlen];
		plink2::PgenHeaderCtrl header_ctrl2;
		uintptr_t pgfi_alloc_ct2 = 0;

		err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr, bind_data->raw_variant_ct,
		                             bind_data->raw_sample_ct, &header_ctrl2, &tmp_pgfi, &pgfi_alloc_ct2, errstr_buf2);
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&tmp_pgfi, &ce);
			throw IOException("read_pfile: failed to init PgenReader for sample-orient pre-read: %s", errstr_buf2);
		}

		AlignedBuffer pgfi_alloc2;
		if (pgfi_alloc_ct2 > 0) {
			pgfi_alloc2.Allocate(pgfi_alloc_ct2 * plink2::kCacheline);
		}

		uint32_t max_vrec2 = 0;
		uintptr_t pgr_alloc_ct2 = 0;
		err = plink2::PgfiInitPhase2(header_ctrl2, 0, 0, 0, 0, tmp_pgfi.raw_variant_ct, &max_vrec2, &tmp_pgfi,
		                             pgfi_alloc2.As<unsigned char>(), &pgr_alloc_ct2, errstr_buf2);
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&tmp_pgfi, &ce);
			throw IOException("read_pfile: failed to init PgenReader for sample-orient pre-read (phase 2): %s",
			                  errstr_buf2);
		}

		plink2::PgenReader tmp_pgr;
		plink2::PreinitPgr(&tmp_pgr);
		AlignedBuffer pgr_alloc2;
		if (pgr_alloc_ct2 > 0) {
			pgr_alloc2.Allocate(pgr_alloc_ct2 * plink2::kCacheline);
		}

		err = plink2::PgrInit(bind_data->pgen_path.c_str(), max_vrec2, &tmp_pgfi, &tmp_pgr,
		                      pgr_alloc2.As<unsigned char>());
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&tmp_pgr, &ce);
			ce = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&tmp_pgfi, &ce);
			throw IOException("read_pfile: PgrInit failed for sample-orient pre-read");
		}

		// Allocate genovec buffer
		uintptr_t genovec_wc = plink2::NypCtToAlignedWordCt(bind_data->raw_sample_ct);
		AlignedBuffer genovec_buf2;
		genovec_buf2.Allocate(genovec_wc * sizeof(uintptr_t));
		std::memset(genovec_buf2.ptr, 0, genovec_wc * sizeof(uintptr_t));

		vector<int8_t> tmp_bytes(bind_data->raw_sample_ct);

		// Set up sample subsetting for pre-read (reuse BuildSampleSubset for consistency)
		SampleSubset preread_subset;
		plink2::PgrSampleSubsetIndex pssi2;
		if (bind_data->has_sample_subset) {
			preread_subset = BuildSampleSubset(bind_data->raw_sample_ct, bind_data->sample_indices);
			plink2::PgrSetSampleSubsetIndex(preread_subset.CumulativePopcounts(), &tmp_pgr, &pssi2);
		} else {
			plink2::PgrClearSampleSubsetIndex(&tmp_pgr, &pssi2);
		}

		// Allocate phase buffers if needed
		AlignedBuffer preread_phasepresent;
		AlignedBuffer preread_phaseinfo;
		vector<int8_t> preread_phased_pairs;
		if (bind_data->include_phased) {
			uintptr_t phase_wc = plink2::BitCtToAlignedWordCt(output_sample_ct);
			preread_phasepresent.Allocate(phase_wc * sizeof(uintptr_t));
			preread_phaseinfo.Allocate(phase_wc * sizeof(uintptr_t));
			preread_phased_pairs.resize(static_cast<size_t>(output_sample_ct) * 2);
		}

		// Allocate dosage buffers if needed
		AlignedBuffer preread_dosage_present;
		AlignedBuffer preread_dosage_main;
		vector<double> preread_dosage_doubles;
		if (bind_data->include_dosages) {
			uintptr_t dosage_present_wc = plink2::BitCtToAlignedWordCt(bind_data->raw_sample_ct);
			preread_dosage_present.Allocate(dosage_present_wc * sizeof(uintptr_t));
			std::memset(preread_dosage_present.ptr, 0, dosage_present_wc * sizeof(uintptr_t));
			preread_dosage_main.Allocate(bind_data->raw_sample_ct * sizeof(uint16_t));
			std::memset(preread_dosage_main.ptr, 0, bind_data->raw_sample_ct * sizeof(uint16_t));
			preread_dosage_doubles.resize(output_sample_ct, 0.0);
		}

		// Pre-read genotypes for each effective variant
		if (bind_data->include_dosages) {
			bind_data->dosage_matrix.resize(effective_variant_ct);
		} else {
			bind_data->genotype_matrix.resize(effective_variant_ct);
		}
		for (uint32_t ev = 0; ev < effective_variant_ct; ev++) {
			uint32_t vidx = bind_data->has_effective_variant_list ? bind_data->effective_variant_indices[ev] : ev;

			const uintptr_t *si_ptr = bind_data->has_sample_subset ? preread_subset.SampleInclude() : nullptr;

			if (bind_data->include_dosages) {
				uint32_t dosage_ct = 0;
				err = plink2::PgrGetD(si_ptr, pssi2, output_sample_ct, vidx, &tmp_pgr, genovec_buf2.As<uintptr_t>(),
				                      preread_dosage_present.As<uintptr_t>(), preread_dosage_main.As<uint16_t>(),
				                      &dosage_ct);
				if (err != plink2::kPglRetSuccess) {
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&tmp_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&tmp_pgfi, &ce);
					throw IOException("read_pfile: PgrGetD failed for variant %u during sample-orient pre-read", vidx);
				}
				plink2::Dosage16ToDoublesMinus9(genovec_buf2.As<uintptr_t>(), preread_dosage_present.As<uintptr_t>(),
				                                preread_dosage_main.As<uint16_t>(), output_sample_ct, dosage_ct,
				                                preread_dosage_doubles.data());
				bind_data->dosage_matrix[ev].assign(preread_dosage_doubles.begin(),
				                                    preread_dosage_doubles.begin() + output_sample_ct);
			} else if (bind_data->include_phased) {
				uint32_t phasepresent_ct = 0;
				err = plink2::PgrGetP(si_ptr, pssi2, output_sample_ct, vidx, &tmp_pgr, genovec_buf2.As<uintptr_t>(),
				                      preread_phasepresent.As<uintptr_t>(), preread_phaseinfo.As<uintptr_t>(),
				                      &phasepresent_ct);
				if (err != plink2::kPglRetSuccess) {
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&tmp_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&tmp_pgfi, &ce);
					throw IOException("read_pfile: PgrGetP failed for variant %u during sample-orient pre-read", vidx);
				}
				plink2::GenoarrToBytesMinus9(genovec_buf2.As<uintptr_t>(), output_sample_ct, tmp_bytes.data());
				UnpackPhasedGenotypes(tmp_bytes.data(), preread_phasepresent.As<uintptr_t>(),
				                      preread_phaseinfo.As<uintptr_t>(), output_sample_ct, preread_phased_pairs.data());
				// Apply genotype_range per-element filtering (phased: check diploid sum)
				if (bind_data->genotype_filter.active && !genotype_range_all_pass.empty() &&
				    !genotype_range_all_pass[ev]) {
					for (uint32_t s = 0; s < output_sample_ct; s++) {
						int8_t a1 = preread_phased_pairs[s * 2];
						int8_t a2 = preread_phased_pairs[s * 2 + 1];
						if (a1 != -9 && !bind_data->genotype_filter.range.Passes(static_cast<double>(a1 + a2))) {
							preread_phased_pairs[s * 2] = -9;
							preread_phased_pairs[s * 2 + 1] = -9;
						}
					}
				}
				bind_data->genotype_matrix[ev].assign(preread_phased_pairs.begin(),
				                                      preread_phased_pairs.begin() + output_sample_ct * 2);
			} else {
				err = plink2::PgrGet(si_ptr, pssi2, output_sample_ct, vidx, &tmp_pgr, genovec_buf2.As<uintptr_t>());
				if (err != plink2::kPglRetSuccess) {
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&tmp_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&tmp_pgfi, &ce);
					throw IOException("read_pfile: PgrGet failed for variant %u during sample-orient pre-read", vidx);
				}
				plink2::GenoarrToBytesMinus9(genovec_buf2.As<uintptr_t>(), output_sample_ct, tmp_bytes.data());
				// Apply genotype_range per-element filtering
				if (bind_data->genotype_filter.active && !genotype_range_all_pass.empty() &&
				    !genotype_range_all_pass[ev]) {
					for (uint32_t s = 0; s < output_sample_ct; s++) {
						int8_t geno = tmp_bytes[s];
						if (geno != -9 && !bind_data->genotype_filter.range.Passes(static_cast<double>(geno))) {
							tmp_bytes[s] = -9;
						}
					}
				}
				bind_data->genotype_matrix[ev].assign(tmp_bytes.begin(), tmp_bytes.begin() + output_sample_ct);
			}
		}

		// Cleanup temporary pgenlib state
		{
			plink2::PglErr ce = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&tmp_pgr, &ce);
			ce = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&tmp_pgfi, &ce);
		}
	} else {
		// Variant mode: same as read_pgen
		uint32_t output_sample_ct = bind_data->OutputSampleCt();

		// Resolve genotype output mode
		string genotypes_str = "auto";
		auto genotypes_it = input.named_parameters.find("genotypes");
		if (genotypes_it != input.named_parameters.end()) {
			genotypes_str = genotypes_it->second.GetValue<string>();
		}
		bind_data->genotype_mode = ResolveGenotypeMode(genotypes_str, output_sample_ct, "read_pfile");

		// Validate incompatible combinations for aggregate modes
		if (IsAggregateGenotypeMode(bind_data->genotype_mode)) {
			if (bind_data->include_phased) {
				throw InvalidInputException("read_pfile: genotypes := '%s' is incompatible with phased := true",
				                            bind_data->genotype_mode == GenotypeMode::COUNTS ? "counts" : "stats");
			}
			if (bind_data->include_dosages) {
				throw InvalidInputException("read_pfile: genotypes := '%s' is incompatible with dosages := true",
				                            bind_data->genotype_mode == GenotypeMode::COUNTS ? "counts" : "stats");
			}
		}

		// Build SampleSubset for COUNTS/STATS mode if needed (for PgrGetCounts in scan)
		if (IsAggregateGenotypeMode(bind_data->genotype_mode) && bind_data->has_sample_subset &&
		    !bind_data->count_filter_subset) {
			bind_data->count_filter_subset =
			    make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, bind_data->sample_indices));
		}

		if (bind_data->genotype_mode == GenotypeMode::COLUMNS) {
			// Columns mode: one scalar column per output sample
			names = {"CHROM", "POS", "ID", "REF", "ALT"};
			return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
			                LogicalType::VARCHAR};

			bind_data->columns_mode_first_geno_col = names.size();
			bind_data->columns_mode_geno_col_count = output_sample_ct;

			LogicalType col_type = bind_data->include_dosages ? LogicalType::DOUBLE : LogicalType::TINYINT;
			for (uint32_t s = 0; s < output_sample_ct; s++) {
				uint32_t file_idx = bind_data->has_sample_subset ? bind_data->sample_indices[s] : s;
				string col_name = bind_data->sample_info.iids[file_idx];
				bind_data->genotype_column_names.push_back(col_name);
				names.push_back(col_name);
				return_types.push_back(col_type);
			}
		} else if (bind_data->genotype_mode == GenotypeMode::STRUCT) {
			// STRUCT mode: single genotypes column with one named field per output sample
			child_list_t<LogicalType> struct_children;
			LogicalType field_type = bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
			                         : bind_data->include_dosages ? LogicalType::DOUBLE
			                                                      : LogicalType::TINYINT;

			for (uint32_t s = 0; s < output_sample_ct; s++) {
				uint32_t file_idx = bind_data->has_sample_subset ? bind_data->sample_indices[s] : s;
				string col_name = bind_data->sample_info.iids[file_idx];
				bind_data->genotype_column_names.push_back(col_name);
				struct_children.push_back({col_name, field_type});
			}

			names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
			return_types = {LogicalType::VARCHAR, LogicalType::INTEGER,
			                LogicalType::VARCHAR, LogicalType::VARCHAR,
			                LogicalType::VARCHAR, LogicalType::STRUCT(std::move(struct_children))};
		} else if (bind_data->genotype_mode == GenotypeMode::COUNTS) {
			names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
			return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
			                LogicalType::VARCHAR, LogicalType::VARCHAR, MakeGenotypeCountsType()};
		} else if (bind_data->genotype_mode == GenotypeMode::STATS) {
			names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
			return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
			                LogicalType::VARCHAR, LogicalType::VARCHAR, MakeGenotypeStatsType()};
		} else {
			LogicalType variant_elem_type = bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
			                                : bind_data->include_dosages ? LogicalType::DOUBLE
			                                                             : LogicalType::TINYINT;
			LogicalType geno_type = bind_data->genotype_mode == GenotypeMode::ARRAY
			                            ? LogicalType::ARRAY(variant_elem_type, output_sample_ct)
			                            : LogicalType::LIST(variant_elem_type);

			names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
			return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
			                LogicalType::VARCHAR, LogicalType::VARCHAR, geno_type};
		}
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
	state->orient_mode = bind_data.orient_mode;
	state->max_threads_config = GetPlinkingMaxThreads(context);

	if (bind_data.orient_mode == OrientMode::GENOTYPE) {
		state->total_count = bind_data.EffectiveVariantCt();
		// Check if genotype column is projected
		state->need_genotypes = false;
		for (auto col_id : input.column_ids) {
			if (col_id == bind_data.geno_genotype_col) {
				state->need_genotypes = true;
				break;
			}
		}
	} else if (bind_data.orient_mode == OrientMode::SAMPLE) {
		state->total_count = bind_data.OutputSampleCt();
		// Check if genotypes column(s) are projected
		state->need_genotypes = false;
		for (auto col_id : input.column_ids) {
			if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}
			if (bind_data.genotype_mode == GenotypeMode::COLUMNS) {
				if (col_id >= bind_data.columns_mode_first_geno_col &&
				    col_id < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
					state->need_genotypes = true;
					break;
				}
			} else if (col_id == bind_data.sample_orient_genotypes_col) {
				state->need_genotypes = true;
				break;
			}
		}
	} else {
		state->total_count = bind_data.EffectiveVariantCt();
		state->need_genotypes = false;
		for (auto col_id : input.column_ids) {
			if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}
			if (bind_data.genotype_mode == GenotypeMode::COLUMNS) {
				if (col_id >= bind_data.columns_mode_first_geno_col &&
				    col_id < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
					state->need_genotypes = true;
					break;
				}
			} else if (col_id == PfileBindData::GENOTYPES_COL) {
				// For COUNTS/STATS we use PgrGetCounts, not PgrGet — don't set need_genotypes
				// (need_genotypes controls whether PgrGet is called in scan)
				if (!IsAggregateGenotypeMode(bind_data.genotype_mode)) {
					state->need_genotypes = true;
				}
				break;
			}
		}
	}

	// need_pgen_reader: true if we need PgrGet (need_genotypes), PgrGetCounts for filters,
	// or PgrGetCounts for COUNTS/STATS aggregate modes in variant orient
	bool need_aggregate_pgen = false;
	if (IsAggregateGenotypeMode(bind_data.genotype_mode) && bind_data.orient_mode == OrientMode::VARIANT) {
		for (auto col_id : input.column_ids) {
			if (col_id == PfileBindData::GENOTYPES_COL) {
				need_aggregate_pgen = true;
				break;
			}
		}
	}
	state->need_pgen_reader = state->need_genotypes || bind_data.count_filter.HasFilter() ||
	                          bind_data.genotype_filter.active || need_aggregate_pgen;

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

	if (!gstate.need_pgen_reader || bind_data.orient_mode == OrientMode::SAMPLE) {
		// Sample-orient mode uses pre-read genotype matrix — no per-thread PgenReader needed
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

	// Allocate phase buffers if phased output requested
	if (bind_data.include_phased) {
		uint32_t output_sample_ct = bind_data.OutputSampleCt();
		uintptr_t phase_wc = plink2::BitCtToAlignedWordCt(output_sample_ct);
		state->phasepresent_buf.Allocate(phase_wc * sizeof(uintptr_t));
		state->phaseinfo_buf.Allocate(phase_wc * sizeof(uintptr_t));
		state->phased_pairs.resize(static_cast<size_t>(output_sample_ct) * 2);
	}

	// Allocate dosage buffers if dosage output requested
	if (bind_data.include_dosages) {
		uint32_t output_sample_ct = bind_data.OutputSampleCt();
		uintptr_t dosage_present_wc = plink2::BitCtToAlignedWordCt(genovec_sample_ct);
		state->dosage_present_buf.Allocate(dosage_present_wc * sizeof(uintptr_t));
		std::memset(state->dosage_present_buf.ptr, 0, dosage_present_wc * sizeof(uintptr_t));

		state->dosage_main_buf.Allocate(genovec_sample_ct * sizeof(uint16_t));
		std::memset(state->dosage_main_buf.ptr, 0, genovec_sample_ct * sizeof(uint16_t));

		state->dosage_doubles.resize(output_sample_ct, 0.0);
	}

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
static constexpr uint32_t PFILE_GENOTYPE_BATCH_SIZE = 64;

static void PfileDefaultScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();
	auto &gstate = data_p.global_state->Cast<PfileGlobalState>();
	auto &lstate = data_p.local_state->Cast<PfileLocalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t total_variants = gstate.total_count;

	uint32_t output_sample_ct = bind_data.OutputSampleCt();

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(PFILE_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_idx.fetch_add(claim_size);
		if (batch_start >= total_variants) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, total_variants);

		for (uint32_t effective_pos = batch_start; effective_pos < batch_end; effective_pos++) {
			uint32_t vidx = ResolveVariantIdx(bind_data, effective_pos);

			// Count filter + genotype range pre-decompression check
			bool geno_range_all_pass = true;
			if ((bind_data.count_filter.HasFilter() || bind_data.genotype_filter.active) && lstate.initialized) {
				STD_ARRAY_DECL(uint32_t, 4, genocounts);
				const uintptr_t *cf_si = (bind_data.has_sample_subset && bind_data.count_filter_subset)
				                             ? bind_data.count_filter_subset->SampleInclude()
				                             : nullptr;
				const uintptr_t *cf_iv = (bind_data.has_sample_subset && bind_data.count_filter_subset)
				                             ? bind_data.count_filter_subset->InterleavedVec()
				                             : nullptr;
				uint32_t cf_sc = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.raw_sample_ct;
				plink2::PglErr cf_err =
				    plink2::PgrGetCounts(cf_si, cf_iv, lstate.pssi, cf_sc, vidx, &lstate.pgr, genocounts);
				if (cf_err != plink2::kPglRetSuccess) {
					throw IOException("read_pfile: PgrGetCounts failed for variant %u", vidx);
				}
				auto pf = CheckPreDecompFilters(bind_data.count_filter, bind_data.genotype_filter, genocounts, cf_sc);
				if (pf.skip) {
					continue;
				}
				geno_range_all_pass = pf.all_pass;
			}

			// Read genotype data if needed
			bool genotypes_read = false;
			if (gstate.need_genotypes && lstate.initialized) {
				const uintptr_t *sample_include =
				    bind_data.has_sample_subset ? lstate.sample_include_buf.As<uintptr_t>() : nullptr;

				if (bind_data.include_dosages) {
					uint32_t dosage_ct = 0;
					plink2::PglErr err =
					    plink2::PgrGetD(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
					                    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
					                    lstate.dosage_main_buf.As<uint16_t>(), &dosage_ct);
					if (err != plink2::kPglRetSuccess) {
						throw IOException("read_pfile: PgrGetD failed for variant %u", vidx);
					}
					plink2::Dosage16ToDoublesMinus9(lstate.genovec_buf.As<uintptr_t>(),
					                                lstate.dosage_present_buf.As<uintptr_t>(),
					                                lstate.dosage_main_buf.As<uint16_t>(), output_sample_ct, dosage_ct,
					                                lstate.dosage_doubles.data());
				} else if (bind_data.include_phased) {
					plink2::PglErr err =
					    plink2::PgrGetP(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
					                    lstate.genovec_buf.As<uintptr_t>(), lstate.phasepresent_buf.As<uintptr_t>(),
					                    lstate.phaseinfo_buf.As<uintptr_t>(), &lstate.phasepresent_ct);
					if (err != plink2::kPglRetSuccess) {
						throw IOException("read_pfile: PgrGetP failed for variant %u", vidx);
					}
					plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
					                             lstate.genotype_bytes.data());
					UnpackPhasedGenotypes(lstate.genotype_bytes.data(), lstate.phasepresent_buf.As<uintptr_t>(),
					                      lstate.phaseinfo_buf.As<uintptr_t>(), output_sample_ct,
					                      lstate.phased_pairs.data());
				} else {
					plink2::PglErr err = plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx,
					                                    &lstate.pgr, lstate.genovec_buf.As<uintptr_t>());
					if (err != plink2::kPglRetSuccess) {
						throw IOException("read_pfile: PgrGet failed for variant %u", vidx);
					}
					plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
					                             lstate.genotype_bytes.data());
				}
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
					auto val = bind_data.variants.GetChrom(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case PfileBindData::POS_COL: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.GetPos(vidx);
					break;
				}
				case PfileBindData::ID_COL: {
					auto val = bind_data.variants.GetId(vidx);
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case PfileBindData::REF_COL: {
					auto val = bind_data.variants.GetRef(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case PfileBindData::ALT_COL: {
					auto val = bind_data.variants.GetAlt(vidx);
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case PfileBindData::GENOTYPES_COL: {
					if (bind_data.genotype_mode == GenotypeMode::COLUMNS) {
						// In COLUMNS mode, file_col 5 is the first genotype column (sample_pos 0)
						if (genotypes_read) {
							idx_t sample_pos = file_col - bind_data.columns_mode_first_geno_col;
							if (bind_data.include_dosages) {
								double dosage = lstate.dosage_doubles[sample_pos];
								if (dosage == -9.0) {
									FlatVector::SetNull(vec, rows_emitted, true);
								} else {
									FlatVector::GetData<double>(vec)[rows_emitted] = dosage;
								}
							} else {
								int8_t geno = lstate.genotype_bytes[sample_pos];
								if (geno == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(geno)))) {
									FlatVector::SetNull(vec, rows_emitted, true);
								} else {
									FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
								}
							}
						} else {
							FlatVector::SetNull(vec, rows_emitted, true);
						}
						break;
					}
					if (bind_data.genotype_mode == GenotypeMode::STRUCT) {
						if (!genotypes_read) {
							FlatVector::SetNull(vec, rows_emitted, true);
							break;
						}
						auto &entries = StructVector::GetEntries(vec);
						for (idx_t s = 0; s < output_sample_ct; s++) {
							auto &child_vec = *entries[s];
							if (bind_data.include_phased) {
								auto &allele_vec = ArrayVector::GetEntry(child_vec);
								auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
								auto &pair_validity = FlatVector::Validity(child_vec);
								idx_t allele_base = rows_emitted * 2;
								int8_t a1 = lstate.phased_pairs[s * 2];
								int8_t a2 = lstate.phased_pairs[s * 2 + 1];
								if (a1 == -9) {
									pair_validity.SetInvalid(rows_emitted);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = a2;
								}
							} else if (bind_data.include_dosages) {
								double dosage = lstate.dosage_doubles[s];
								if (dosage == -9.0) {
									FlatVector::SetNull(child_vec, rows_emitted, true);
								} else {
									FlatVector::GetData<double>(child_vec)[rows_emitted] = dosage;
								}
							} else {
								int8_t geno = lstate.genotype_bytes[s];
								if (geno == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(geno)))) {
									FlatVector::SetNull(child_vec, rows_emitted, true);
								} else {
									FlatVector::GetData<int8_t>(child_vec)[rows_emitted] = geno;
								}
							}
						}
						break;
					}
					if (IsAggregateGenotypeMode(bind_data.genotype_mode)) {
						// COUNTS/STATS: use PgrGetCounts (no decompression needed)
						if (!lstate.initialized) {
							FlatVector::SetNull(vec, rows_emitted, true);
							break;
						}
						STD_ARRAY_DECL(uint32_t, 4, genocounts);
						const uintptr_t *agg_si = (bind_data.has_sample_subset && bind_data.count_filter_subset)
						                              ? bind_data.count_filter_subset->SampleInclude()
						                              : nullptr;
						const uintptr_t *agg_iv = (bind_data.has_sample_subset && bind_data.count_filter_subset)
						                              ? bind_data.count_filter_subset->InterleavedVec()
						                              : nullptr;
						uint32_t agg_sc =
						    bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.raw_sample_ct;

						plink2::PglErr agg_err =
						    plink2::PgrGetCounts(agg_si, agg_iv, lstate.pssi, agg_sc, vidx, &lstate.pgr, genocounts);
						if (agg_err != plink2::kPglRetSuccess) {
							throw IOException("read_pfile: PgrGetCounts failed for variant %u", vidx);
						}

						auto &entries = StructVector::GetEntries(vec);
						FlatVector::GetData<uint32_t>(*entries[0])[rows_emitted] = genocounts[0];
						FlatVector::GetData<uint32_t>(*entries[1])[rows_emitted] = genocounts[1];
						FlatVector::GetData<uint32_t>(*entries[2])[rows_emitted] = genocounts[2];
						FlatVector::GetData<uint32_t>(*entries[3])[rows_emitted] = genocounts[3];

						if (bind_data.genotype_mode == GenotypeMode::STATS) {
							uint32_t n = genocounts[0] + genocounts[1] + genocounts[2];
							uint32_t total = n + genocounts[3];
							FlatVector::GetData<uint32_t>(*entries[4])[rows_emitted] = n;
							if (n == 0) {
								FlatVector::GetData<double>(*entries[5])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
								FlatVector::GetData<double>(*entries[6])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
								FlatVector::GetData<double>(*entries[9])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
							} else {
								double af = (static_cast<double>(genocounts[1]) + 2.0 * genocounts[2]) / (2.0 * n);
								FlatVector::GetData<double>(*entries[5])[rows_emitted] = af;
								FlatVector::GetData<double>(*entries[6])[rows_emitted] = std::min(af, 1.0 - af);
								FlatVector::GetData<double>(*entries[9])[rows_emitted] =
								    static_cast<double>(genocounts[1]) / static_cast<double>(n);
							}
							if (total == 0) {
								FlatVector::GetData<double>(*entries[7])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
							} else {
								FlatVector::GetData<double>(*entries[7])[rows_emitted] =
								    static_cast<double>(genocounts[3]) / static_cast<double>(total);
							}
							FlatVector::GetData<uint32_t>(*entries[8])[rows_emitted] = genocounts[1] + genocounts[2];
						}
						break;
					}
					if (!genotypes_read) {
						if (bind_data.genotype_mode == GenotypeMode::LIST) {
							auto *list_data = FlatVector::GetData<list_entry_t>(vec);
							list_data[rows_emitted].offset = ListVector::GetListSize(vec);
							list_data[rows_emitted].length = 0;
						}
						FlatVector::SetNull(vec, rows_emitted, true);
						break;
					}

					if (bind_data.include_phased) {
						if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
							auto &pair_vec = ArrayVector::GetEntry(vec);
							auto &allele_vec = ArrayVector::GetEntry(pair_vec);
							auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
							auto &pair_validity = FlatVector::Validity(pair_vec);

							idx_t pair_base = rows_emitted * static_cast<idx_t>(output_sample_ct);
							for (idx_t s = 0; s < output_sample_ct; s++) {
								idx_t pair_idx = pair_base + s;
								idx_t allele_base = pair_idx * 2;
								int8_t a1 = lstate.phased_pairs[s * 2];
								int8_t a2 = lstate.phased_pairs[s * 2 + 1];
								if (a1 == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(a1 + a2)))) {
									pair_validity.SetInvalid(pair_idx);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = a2;
								}
							}
						} else {
							auto list_offset = ListVector::GetListSize(vec);
							ListVector::Reserve(vec, list_offset + output_sample_ct);
							auto &pair_vec = ListVector::GetEntry(vec);
							auto &allele_vec = ArrayVector::GetEntry(pair_vec);
							auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
							auto &pair_validity = FlatVector::Validity(pair_vec);

							for (idx_t s = 0; s < output_sample_ct; s++) {
								idx_t pair_idx = list_offset + s;
								idx_t allele_base = pair_idx * 2;
								int8_t a1 = lstate.phased_pairs[s * 2];
								int8_t a2 = lstate.phased_pairs[s * 2 + 1];
								if (a1 == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(a1 + a2)))) {
									pair_validity.SetInvalid(pair_idx);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = a2;
								}
							}

							auto *list_data = FlatVector::GetData<list_entry_t>(vec);
							list_data[rows_emitted].offset = list_offset;
							list_data[rows_emitted].length = output_sample_ct;
							ListVector::SetListSize(vec, list_offset + output_sample_ct);
						}
					} else if (bind_data.include_dosages) {
						// Dosage output: ARRAY(DOUBLE, N) or LIST(DOUBLE)
						if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
							auto array_size = static_cast<idx_t>(output_sample_ct);
							auto &child = ArrayVector::GetEntry(vec);
							auto *child_data = FlatVector::GetData<double>(child);
							auto &child_validity = FlatVector::Validity(child);

							idx_t base = rows_emitted * array_size;
							for (idx_t s = 0; s < array_size; s++) {
								double dosage = lstate.dosage_doubles[s];
								if (dosage == -9.0) {
									child_validity.SetInvalid(base + s);
									child_data[base + s] = 0.0;
								} else {
									child_data[base + s] = dosage;
								}
							}
						} else {
							// LIST(DOUBLE)
							auto list_offset = ListVector::GetListSize(vec);
							ListVector::Reserve(vec, list_offset + output_sample_ct);
							auto &child = ListVector::GetEntry(vec);
							auto *child_data = FlatVector::GetData<double>(child);
							auto &child_validity = FlatVector::Validity(child);
							for (idx_t s = 0; s < output_sample_ct; s++) {
								double dosage = lstate.dosage_doubles[s];
								if (dosage == -9.0) {
									child_validity.SetInvalid(list_offset + s);
									child_data[list_offset + s] = 0.0;
								} else {
									child_data[list_offset + s] = dosage;
								}
							}
							auto *list_data = FlatVector::GetData<list_entry_t>(vec);
							list_data[rows_emitted].offset = list_offset;
							list_data[rows_emitted].length = output_sample_ct;
							ListVector::SetListSize(vec, list_offset + output_sample_ct);
						}
					} else {
						if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
							auto array_size = static_cast<idx_t>(output_sample_ct);
							auto &child = ArrayVector::GetEntry(vec);
							auto *child_data = FlatVector::GetData<int8_t>(child);
							auto &child_validity = FlatVector::Validity(child);

							idx_t base = rows_emitted * array_size;
							for (idx_t s = 0; s < array_size; s++) {
								int8_t geno = lstate.genotype_bytes[s];
								if (geno == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(geno)))) {
									child_validity.SetInvalid(base + s);
									child_data[base + s] = 0;
								} else {
									child_data[base + s] = geno;
								}
							}
						} else {
							auto list_offset = ListVector::GetListSize(vec);
							ListVector::Reserve(vec, list_offset + output_sample_ct);
							auto &child = ListVector::GetEntry(vec);
							auto *child_data = FlatVector::GetData<int8_t>(child);
							auto &child_validity = FlatVector::Validity(child);
							for (idx_t s = 0; s < output_sample_ct; s++) {
								int8_t geno = lstate.genotype_bytes[s];
								if (geno == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(geno)))) {
									child_validity.SetInvalid(list_offset + s);
									child_data[list_offset + s] = 0;
								} else {
									child_data[list_offset + s] = geno;
								}
							}
							auto *list_data = FlatVector::GetData<list_entry_t>(vec);
							list_data[rows_emitted].offset = list_offset;
							list_data[rows_emitted].length = output_sample_ct;
							ListVector::SetListSize(vec, list_offset + output_sample_ct);
						}
					}
					break;
				}
				default: {
					// Columns mode: individual genotype columns
					if (bind_data.genotype_mode == GenotypeMode::COLUMNS &&
					    file_col >= bind_data.columns_mode_first_geno_col &&
					    file_col < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
						if (genotypes_read) {
							idx_t sample_pos = file_col - bind_data.columns_mode_first_geno_col;
							if (bind_data.include_dosages) {
								double dosage = lstate.dosage_doubles[sample_pos];
								if (dosage == -9.0) {
									FlatVector::SetNull(vec, rows_emitted, true);
								} else {
									FlatVector::GetData<double>(vec)[rows_emitted] = dosage;
								}
							} else {
								int8_t geno = lstate.genotype_bytes[sample_pos];
								if (geno == -9 ||
								    (bind_data.genotype_filter.active && !geno_range_all_pass &&
								     !bind_data.genotype_filter.range.Passes(static_cast<double>(geno)))) {
									FlatVector::SetNull(vec, rows_emitted, true);
								} else {
									FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
								}
							}
						} else {
							FlatVector::SetNull(vec, rows_emitted, true);
						}
					}
					break;
				}
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

	auto &column_ids = gstate.column_ids;
	uint32_t total_effective_variants = gstate.total_count;
	uint32_t output_sample_ct = bind_data.OutputSampleCt();

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		// Claim a new batch if current batch is exhausted
		if (lstate.batch_exhausted) {
			lstate.batch_start = gstate.next_idx.fetch_add(PFILE_GENOTYPE_BATCH_SIZE);
			if (lstate.batch_start >= total_effective_variants) {
				break; // no more work
			}
			lstate.batch_end = std::min(lstate.batch_start + PFILE_GENOTYPE_BATCH_SIZE, total_effective_variants);
			lstate.current_variant_in_batch = lstate.batch_start;
			lstate.current_sample_in_variant = 0;
			lstate.batch_variant_loaded = false;
			lstate.batch_exhausted = false;
		}

		// Check if we've exhausted the current batch
		if (lstate.current_variant_in_batch >= lstate.batch_end) {
			lstate.batch_exhausted = true;
			continue;
		}

		uint32_t vidx = ResolveVariantIdx(bind_data, lstate.current_variant_in_batch);

		// Count filter + genotype range pre-decompression check
		if ((bind_data.count_filter.HasFilter() || bind_data.genotype_filter.active) && lstate.initialized &&
		    !lstate.batch_variant_loaded) {
			STD_ARRAY_DECL(uint32_t, 4, genocounts);
			const uintptr_t *cf_si = (bind_data.has_sample_subset && bind_data.count_filter_subset)
			                             ? bind_data.count_filter_subset->SampleInclude()
			                             : nullptr;
			const uintptr_t *cf_iv = (bind_data.has_sample_subset && bind_data.count_filter_subset)
			                             ? bind_data.count_filter_subset->InterleavedVec()
			                             : nullptr;
			uint32_t cf_sc = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.raw_sample_ct;
			plink2::PglErr cf_err =
			    plink2::PgrGetCounts(cf_si, cf_iv, lstate.pssi, cf_sc, vidx, &lstate.pgr, genocounts);
			if (cf_err != plink2::kPglRetSuccess) {
				throw IOException("read_pfile: PgrGetCounts failed for variant %u", vidx);
			}
			auto pf = CheckPreDecompFilters(bind_data.count_filter, bind_data.genotype_filter, genocounts, cf_sc);
			if (pf.skip) {
				lstate.current_variant_in_batch++;
				lstate.current_sample_in_variant = 0;
				continue;
			}
			lstate.geno_range_all_pass = pf.all_pass;
		}

		// Load genotypes for current variant if not yet loaded
		if (!lstate.batch_variant_loaded && gstate.need_genotypes && lstate.initialized) {
			const uintptr_t *sample_include =
			    bind_data.has_sample_subset ? lstate.sample_include_buf.As<uintptr_t>() : nullptr;

			if (bind_data.include_dosages) {
				uint32_t dosage_ct = 0;
				plink2::PglErr err =
				    plink2::PgrGetD(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
				                    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
				                    lstate.dosage_main_buf.As<uint16_t>(), &dosage_ct);
				if (err != plink2::kPglRetSuccess) {
					throw IOException("read_pfile: PgrGetD failed for variant %u", vidx);
				}
				plink2::Dosage16ToDoublesMinus9(
				    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
				    lstate.dosage_main_buf.As<uint16_t>(), output_sample_ct, dosage_ct, lstate.dosage_doubles.data());
			} else if (bind_data.include_phased) {
				plink2::PglErr err =
				    plink2::PgrGetP(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
				                    lstate.genovec_buf.As<uintptr_t>(), lstate.phasepresent_buf.As<uintptr_t>(),
				                    lstate.phaseinfo_buf.As<uintptr_t>(), &lstate.phasepresent_ct);
				if (err != plink2::kPglRetSuccess) {
					throw IOException("read_pfile: PgrGetP failed for variant %u", vidx);
				}
				plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
				                             lstate.genotype_bytes.data());
				UnpackPhasedGenotypes(lstate.genotype_bytes.data(), lstate.phasepresent_buf.As<uintptr_t>(),
				                      lstate.phaseinfo_buf.As<uintptr_t>(), output_sample_ct,
				                      lstate.phased_pairs.data());
			} else {
				plink2::PglErr err = plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
				                                    lstate.genovec_buf.As<uintptr_t>());
				if (err != plink2::kPglRetSuccess) {
					throw IOException("read_pfile: PgrGet failed for variant %u", vidx);
				}
				plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
				                             lstate.genotype_bytes.data());
			}
			lstate.batch_variant_loaded = true;
		}

		// Emit rows for samples within current variant
		while (lstate.current_sample_in_variant < output_sample_ct && rows_emitted < STANDARD_VECTOR_SIZE) {
			// When genotype_range is active, skip rows for missing/out-of-range genotypes
			if (bind_data.genotype_filter.active && gstate.need_genotypes && lstate.batch_variant_loaded) {
				if (bind_data.include_phased) {
					int8_t a1 = lstate.phased_pairs[lstate.current_sample_in_variant * 2];
					if (a1 == -9) {
						lstate.current_sample_in_variant++;
						continue;
					}
					if (!lstate.geno_range_all_pass) {
						int8_t a2 = lstate.phased_pairs[lstate.current_sample_in_variant * 2 + 1];
						if (!bind_data.genotype_filter.range.Passes(static_cast<double>(a1 + a2))) {
							lstate.current_sample_in_variant++;
							continue;
						}
					}
				} else if (!bind_data.include_dosages) {
					int8_t geno = lstate.genotype_bytes[lstate.current_sample_in_variant];
					if (geno == -9) {
						lstate.current_sample_in_variant++;
						continue;
					}
					if (!lstate.geno_range_all_pass &&
					    !bind_data.genotype_filter.range.Passes(static_cast<double>(geno))) {
						lstate.current_sample_in_variant++;
						continue;
					}
				}
			}

			// Map scan-order sample index to file-order sample index.
			// When subsetting, sample_indices is sorted to match pgenlib's
			// ascending-bit-order output, so genotype_bytes[i] corresponds
			// to sample_indices[i].
			uint32_t sample_file_idx = bind_data.has_sample_subset
			                               ? bind_data.sample_indices[lstate.current_sample_in_variant]
			                               : lstate.current_sample_in_variant;

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
						auto val = bind_data.variants.GetChrom(vidx);
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						break;
					}
					case PfileBindData::POS_COL: {
						FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.GetPos(vidx);
						break;
					}
					case PfileBindData::ID_COL: {
						auto val = bind_data.variants.GetId(vidx);
						if (val.empty()) {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						}
						break;
					}
					case PfileBindData::REF_COL: {
						auto val = bind_data.variants.GetRef(vidx);
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						break;
					}
					case PfileBindData::ALT_COL: {
						auto val = bind_data.variants.GetAlt(vidx);
						if (val.empty() || val == ".") {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
						}
						break;
					}
					}
				} else if (file_col == bind_data.geno_genotype_col) {
					// Genotype column
					if (gstate.need_genotypes && lstate.batch_variant_loaded) {
						if (bind_data.include_phased) {
							// ARRAY(TINYINT, 2) per row
							auto &allele_vec = ArrayVector::GetEntry(vec);
							auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
							idx_t allele_base = rows_emitted * 2;
							int8_t a1 = lstate.phased_pairs[lstate.current_sample_in_variant * 2];
							if (a1 == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
								allele_data[allele_base] = 0;
								allele_data[allele_base + 1] = 0;
							} else {
								allele_data[allele_base] = a1;
								allele_data[allele_base + 1] =
								    lstate.phased_pairs[lstate.current_sample_in_variant * 2 + 1];
							}
						} else if (bind_data.include_dosages) {
							// Scalar DOUBLE
							double dosage = lstate.dosage_doubles[lstate.current_sample_in_variant];
							if (dosage == -9.0) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<double>(vec)[rows_emitted] = dosage;
							}
						} else {
							// Scalar TINYINT
							int8_t geno = lstate.genotype_bytes[lstate.current_sample_in_variant];
							if (geno == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
							}
						}
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
				} else if (file_col >= bind_data.geno_sample_col_start && file_col < bind_data.geno_genotype_col) {
					// Sample metadata column
					idx_t sample_col_rel = file_col - bind_data.geno_sample_col_start;
					idx_t psam_col_idx = bind_data.geno_sample_col_to_psam_col[sample_col_rel];

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
			lstate.current_sample_in_variant++;
		}

		// Advance to next variant if all samples emitted
		if (lstate.current_sample_in_variant >= output_sample_ct) {
			lstate.current_variant_in_batch++;
			lstate.current_sample_in_variant = 0;
			lstate.batch_variant_loaded = false;
		}
	}

	output.SetCardinality(rows_emitted);
}

// ---------------------------------------------------------------------------
// Scan: Sample-orient mode (one row per sample)
// ---------------------------------------------------------------------------

static void PfileSampleOrientScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();
	auto &gstate = data_p.global_state->Cast<PfileGlobalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t total_samples = gstate.total_count;
	uint32_t effective_variant_ct = bind_data.EffectiveVariantCt();

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(PFILE_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_idx.fetch_add(claim_size);
		if (batch_start >= total_samples) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, total_samples);

		for (uint32_t sample_pos = batch_start; sample_pos < batch_end; sample_pos++) {
			// Map scan-order sample index to file-order sample index
			uint32_t sample_file_idx = bind_data.has_sample_subset ? bind_data.sample_indices[sample_pos] : sample_pos;

			for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
				auto file_col = column_ids[out_col];
				if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
					continue;
				}

				auto &vec = output.data[out_col];

				if (bind_data.genotype_mode == GenotypeMode::COLUMNS &&
				    file_col >= bind_data.columns_mode_first_geno_col &&
				    file_col < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
					// Columns mode: individual variant genotype columns
					if (gstate.need_genotypes) {
						idx_t variant_pos = file_col - bind_data.columns_mode_first_geno_col;
						if (bind_data.include_dosages) {
							double dosage = bind_data.dosage_matrix[variant_pos][sample_pos];
							if (dosage == -9.0) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<double>(vec)[rows_emitted] = dosage;
							}
						} else {
							int8_t geno = bind_data.genotype_matrix[variant_pos][sample_pos];
							if (geno == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
							}
						}
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
				} else if (bind_data.genotype_mode != GenotypeMode::COLUMNS &&
				           file_col == bind_data.sample_orient_genotypes_col) {
					// Genotypes column
					if (!gstate.need_genotypes) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else if (bind_data.genotype_mode == GenotypeMode::STRUCT) {
						// STRUCT mode: one field per variant, values from genotype_matrix
						auto &entries = StructVector::GetEntries(vec);
						for (idx_t v = 0; v < effective_variant_ct; v++) {
							auto &child_vec = *entries[v];
							if (bind_data.include_phased) {
								auto &allele_vec = ArrayVector::GetEntry(child_vec);
								auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
								auto &pair_validity = FlatVector::Validity(child_vec);
								idx_t allele_base = rows_emitted * 2;
								int8_t a1 = bind_data.genotype_matrix[v][sample_pos * 2];
								if (a1 == -9) {
									pair_validity.SetInvalid(rows_emitted);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = bind_data.genotype_matrix[v][sample_pos * 2 + 1];
								}
							} else if (bind_data.include_dosages) {
								double dosage = bind_data.dosage_matrix[v][sample_pos];
								if (dosage == -9.0) {
									FlatVector::SetNull(child_vec, rows_emitted, true);
								} else {
									FlatVector::GetData<double>(child_vec)[rows_emitted] = dosage;
								}
							} else {
								int8_t geno = bind_data.genotype_matrix[v][sample_pos];
								if (geno == -9) {
									FlatVector::SetNull(child_vec, rows_emitted, true);
								} else {
									FlatVector::GetData<int8_t>(child_vec)[rows_emitted] = geno;
								}
							}
						}
					} else if (IsAggregateGenotypeMode(bind_data.genotype_mode)) {
						// COUNTS/STATS: accumulate from genotype_matrix across variants
						uint32_t hom_ref = 0, het = 0, hom_alt = 0, missing = 0;
						for (idx_t v = 0; v < effective_variant_ct; v++) {
							int8_t geno = bind_data.genotype_matrix[v][sample_pos];
							switch (geno) {
							case 0:
								hom_ref++;
								break;
							case 1:
								het++;
								break;
							case 2:
								hom_alt++;
								break;
							default:
								missing++;
								break;
							}
						}
						auto &entries = StructVector::GetEntries(vec);
						FlatVector::GetData<uint32_t>(*entries[0])[rows_emitted] = hom_ref;
						FlatVector::GetData<uint32_t>(*entries[1])[rows_emitted] = het;
						FlatVector::GetData<uint32_t>(*entries[2])[rows_emitted] = hom_alt;
						FlatVector::GetData<uint32_t>(*entries[3])[rows_emitted] = missing;

						if (bind_data.genotype_mode == GenotypeMode::STATS) {
							uint32_t n = hom_ref + het + hom_alt;
							uint32_t total = n + missing;
							FlatVector::GetData<uint32_t>(*entries[4])[rows_emitted] = n;
							if (n == 0) {
								FlatVector::GetData<double>(*entries[5])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
								FlatVector::GetData<double>(*entries[6])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
								FlatVector::GetData<double>(*entries[9])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
							} else {
								double af = (static_cast<double>(het) + 2.0 * hom_alt) / (2.0 * n);
								FlatVector::GetData<double>(*entries[5])[rows_emitted] = af;
								FlatVector::GetData<double>(*entries[6])[rows_emitted] = std::min(af, 1.0 - af);
								FlatVector::GetData<double>(*entries[9])[rows_emitted] =
								    static_cast<double>(het) / static_cast<double>(n);
							}
							if (total == 0) {
								FlatVector::GetData<double>(*entries[7])[rows_emitted] =
								    std::numeric_limits<double>::quiet_NaN();
							} else {
								FlatVector::GetData<double>(*entries[7])[rows_emitted] =
								    static_cast<double>(missing) / static_cast<double>(total);
							}
							FlatVector::GetData<uint32_t>(*entries[8])[rows_emitted] = het + hom_alt;
						}
					} else if (bind_data.include_phased) {
						// Phased sample-orient: each variant has stride-2 layout in genotype_matrix
						if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
							auto &pair_vec = ArrayVector::GetEntry(vec);
							auto &allele_vec = ArrayVector::GetEntry(pair_vec);
							auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
							auto &pair_validity = FlatVector::Validity(pair_vec);

							idx_t pair_base = rows_emitted * static_cast<idx_t>(effective_variant_ct);
							for (idx_t v = 0; v < effective_variant_ct; v++) {
								idx_t pair_idx = pair_base + v;
								idx_t allele_base = pair_idx * 2;
								int8_t a1 = bind_data.genotype_matrix[v][sample_pos * 2];
								if (a1 == -9) {
									pair_validity.SetInvalid(pair_idx);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = bind_data.genotype_matrix[v][sample_pos * 2 + 1];
								}
							}
						} else {
							// LIST(ARRAY(TINYINT, 2))
							auto list_offset = ListVector::GetListSize(vec);
							ListVector::Reserve(vec, list_offset + effective_variant_ct);
							auto &pair_vec = ListVector::GetEntry(vec);
							auto &allele_vec = ArrayVector::GetEntry(pair_vec);
							auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
							auto &pair_validity = FlatVector::Validity(pair_vec);

							for (idx_t v = 0; v < effective_variant_ct; v++) {
								idx_t pair_idx = list_offset + v;
								idx_t allele_base = pair_idx * 2;
								int8_t a1 = bind_data.genotype_matrix[v][sample_pos * 2];
								if (a1 == -9) {
									pair_validity.SetInvalid(pair_idx);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = bind_data.genotype_matrix[v][sample_pos * 2 + 1];
								}
							}

							auto *list_data = FlatVector::GetData<list_entry_t>(vec);
							list_data[rows_emitted].offset = list_offset;
							list_data[rows_emitted].length = effective_variant_ct;
							ListVector::SetListSize(vec, list_offset + effective_variant_ct);
						}
					} else if (bind_data.include_dosages) {
						// Dosage sample-orient: ARRAY(DOUBLE, M) or LIST(DOUBLE)
						if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
							auto array_size = static_cast<idx_t>(effective_variant_ct);
							auto &child = ArrayVector::GetEntry(vec);
							auto *child_data = FlatVector::GetData<double>(child);
							auto &child_validity = FlatVector::Validity(child);

							idx_t base = rows_emitted * array_size;
							for (idx_t v = 0; v < array_size; v++) {
								double dosage = bind_data.dosage_matrix[v][sample_pos];
								if (dosage == -9.0) {
									child_validity.SetInvalid(base + v);
									child_data[base + v] = 0.0;
								} else {
									child_data[base + v] = dosage;
								}
							}
						} else {
							// LIST(DOUBLE)
							auto list_offset = ListVector::GetListSize(vec);
							ListVector::Reserve(vec, list_offset + effective_variant_ct);
							auto &child = ListVector::GetEntry(vec);
							auto *child_data = FlatVector::GetData<double>(child);
							auto &child_validity = FlatVector::Validity(child);
							for (idx_t v = 0; v < effective_variant_ct; v++) {
								double dosage = bind_data.dosage_matrix[v][sample_pos];
								if (dosage == -9.0) {
									child_validity.SetInvalid(list_offset + v);
									child_data[list_offset + v] = 0.0;
								} else {
									child_data[list_offset + v] = dosage;
								}
							}
							auto *list_data = FlatVector::GetData<list_entry_t>(vec);
							list_data[rows_emitted].offset = list_offset;
							list_data[rows_emitted].length = effective_variant_ct;
							ListVector::SetListSize(vec, list_offset + effective_variant_ct);
						}
					} else if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
						auto array_size = static_cast<idx_t>(effective_variant_ct);
						auto &child = ArrayVector::GetEntry(vec);
						auto *child_data = FlatVector::GetData<int8_t>(child);
						auto &child_validity = FlatVector::Validity(child);

						idx_t base = rows_emitted * array_size;
						for (idx_t v = 0; v < array_size; v++) {
							int8_t geno = bind_data.genotype_matrix[v][sample_pos];
							if (geno == -9) {
								child_validity.SetInvalid(base + v);
								child_data[base + v] = 0;
							} else {
								child_data[base + v] = geno;
							}
						}
					} else {
						// LIST path
						auto list_offset = ListVector::GetListSize(vec);
						ListVector::Reserve(vec, list_offset + effective_variant_ct);
						auto &child = ListVector::GetEntry(vec);
						auto *child_data = FlatVector::GetData<int8_t>(child);
						auto &child_validity = FlatVector::Validity(child);
						for (idx_t v = 0; v < effective_variant_ct; v++) {
							int8_t geno = bind_data.genotype_matrix[v][sample_pos];
							if (geno == -9) {
								child_validity.SetInvalid(list_offset + v);
								child_data[list_offset + v] = 0;
							} else {
								child_data[list_offset + v] = geno;
							}
						}
						auto *list_data = FlatVector::GetData<list_entry_t>(vec);
						list_data[rows_emitted].offset = list_offset;
						list_data[rows_emitted].length = effective_variant_ct;
						ListVector::SetListSize(vec, list_offset + effective_variant_ct);
					}
				} else {
					// Sample metadata column
					idx_t psam_col_idx = bind_data.sample_orient_col_to_psam_col[file_col];

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
		}
	}

	output.SetCardinality(rows_emitted);
}

// ---------------------------------------------------------------------------
// Scan dispatch
// ---------------------------------------------------------------------------

static void PfileScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();

	switch (bind_data.orient_mode) {
	case OrientMode::GENOTYPE:
		PfileTidyScan(context, data_p, output);
		break;
	case OrientMode::SAMPLE:
		PfileSampleOrientScan(context, data_p, output);
		break;
	default:
		PfileDefaultScan(context, data_p, output);
		break;
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
	read_pfile.named_parameters["orient"] = LogicalType::VARCHAR;
	read_pfile.named_parameters["dosages"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["phased"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["region"] = LogicalType::VARCHAR;
	// Accept ANY for samples and variants — type dispatch handled in bind
	read_pfile.named_parameters["samples"] = LogicalType::ANY;
	read_pfile.named_parameters["variants"] = LogicalType::ANY;
	read_pfile.named_parameters["genotypes"] = LogicalType::VARCHAR;
	read_pfile.named_parameters["af_range"] = LogicalType::ANY;
	read_pfile.named_parameters["ac_range"] = LogicalType::ANY;
	read_pfile.named_parameters["genotype_range"] = LogicalType::ANY;

	loader.RegisterFunction(read_pfile);
}

} // namespace duckdb
