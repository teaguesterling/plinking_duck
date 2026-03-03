#include "pfile_reader.hpp"
#include "plink_common.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"

#include <pgenlib_read.h>
#include <pgenlib_ffi_support.h>

#include <algorithm>
#include <atomic>
#include <cstring>
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
//! Also populates sample_info to avoid a separate file read for LoadSampleInfo.
static PfileSampleMetadata LoadPfileSampleMetadata(ClientContext &context, const string &path,
                                                   SampleInfo &sample_info_out) {
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

		// Extract sample info for IID lookups
		if (iid_idx < fields.size()) {
			const auto &iid = fields[iid_idx];
			if (sample_info_out.iid_to_idx.count(iid)) {
				throw IOException("read_pfile: file '%s' line %llu has duplicate IID '%s'", path,
				                  static_cast<unsigned long long>(i + 1), iid);
			}
			sample_info_out.iids.push_back(iid);
			if (has_fid && fid_idx < fields.size()) {
				sample_info_out.fids.push_back(fields[fid_idx]);
			}
			sample_info_out.iid_to_idx[iid] = sample_info_out.iids.size() - 1;
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

	// Full sample metadata (for tidy mode column output)
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

	// Effective variant list (after region + variant filter intersection)
	// If empty and no filters active, scan all variants sequentially.
	// If non-empty, these are the specific variant indices to scan.
	bool has_effective_variant_list = false;
	vector<uint32_t> effective_variant_indices;

	// Sample-orient mode: pre-read genotype matrix (variant × sample)
	// genotype_matrix[effective_vidx][sample_idx] = genotype value (-9 = missing)
	vector<vector<int8_t>> genotype_matrix;

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
	// For variant/sample-orient mode: atomic counter
	std::atomic<uint32_t> next_idx {0};
	uint32_t total_count = 0; // variants (variant/genotype mode) or samples (sample mode)

	// Projection
	bool need_genotypes = false;
	vector<column_t> column_ids;
	OrientMode orient_mode = OrientMode::VARIANT;

	idx_t MaxThreads() const override {
		// Genotype mode must be single-threaded: the state machine tracks
		// current_variant and current_sample, which are not thread-safe.
		if (orient_mode == OrientMode::GENOTYPE) {
			return 1;
		}
		// Variant and sample modes support parallel scan
		return std::min<idx_t>(total_count / 1000 + 1, 16);
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
	string orient_str;
	bool tidy_flag = false;
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pgen") {
			bind_data->pgen_path = kv.second.GetValue<string>();
		} else if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "tidy") {
			tidy_flag = kv.second.GetValue<bool>();
		} else if (kv.first == "orient") {
			orient_str = kv.second.GetValue<string>();
		} else if (kv.first == "dosages") {
			bind_data->include_dosages = kv.second.GetValue<bool>();
		} else if (kv.first == "phased") {
			bind_data->include_phased = kv.second.GetValue<bool>();
		} else if (kv.first == "region") {
			bind_data->region = ParseRegion(kv.second.GetValue<string>());
		}
		// samples, variants, and genotypes handled after pgenlib init
	}

	bind_data->orient_mode = ResolveOrientMode(orient_str, tidy_flag, "read_pfile");

	if (bind_data->include_dosages) {
		throw NotImplementedException("read_pfile: dosages support is not yet implemented");
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
			bind_data->pvar_path = FindCompanionFile(fs, prefix + ".pgen", {".pvar", ".bim"});
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
			bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
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
			bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
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
	bind_data->variants = LoadVariantMetadataIndex(context, bind_data->pvar_path, "read_pfile");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("read_pfile: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info ---
	// In genotype/sample orient modes, LoadPfileSampleMetadata populates both
	// sample_info and sample_metadata from a single file read. Otherwise, use LoadSampleInfo.
	if (bind_data->orient_mode == OrientMode::GENOTYPE || bind_data->orient_mode == OrientMode::SAMPLE) {
		bind_data->sample_metadata = LoadPfileSampleMetadata(context, bind_data->psam_path, bind_data->sample_info);
	} else {
		bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);
	}

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
			// Build ID→index map from VariantMetadataIndex
			unordered_map<string, uint32_t> id_to_idx;
			for (idx_t i = 0; i < bind_data->variants.variant_ct; i++) {
				auto id = bind_data->variants.GetId(i);
				if (!id.empty()) {
					id_to_idx[id] = static_cast<uint32_t>(i);
				}
			}

			auto &children = ListValue::GetChildren(variants_val);
			for (auto &child : children) {
				auto vid = child.GetValue<string>();
				auto it = id_to_idx.find(vid);
				if (it == id_to_idx.end()) {
					throw InvalidInputException("read_pfile: variant '%s' not found in .pvar", vid);
				}
				bind_data->variant_indices.push_back(it->second);
			}
		} else {
			throw InvalidInputException("read_pfile: variants parameter must be LIST(VARCHAR) or LIST(INTEGER)");
		}

		// Validate no duplicate variant indices (consistent with samples behavior)
		{
			std::unordered_set<uint32_t> seen;
			for (auto idx : bind_data->variant_indices) {
				if (!seen.insert(idx).second) {
					throw InvalidInputException("read_pfile: duplicate variant index %u in variants list", idx);
				}
			}
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
					if (bind_data->variants.GetChrom(vidx) != bind_data->region.chrom) {
						continue;
					}
					int64_t pos = bind_data->variants.GetPos(vidx);
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
	if (bind_data->orient_mode == OrientMode::GENOTYPE) {
		// Check for incompatible genotypes := 'columns' before building schema
		auto genotypes_it_geno = input.named_parameters.find("genotypes");
		if (genotypes_it_geno != input.named_parameters.end()) {
			auto gval = StringUtil::Lower(genotypes_it_geno->second.GetValue<string>());
			if (gval == "columns") {
				throw InvalidInputException(
				    "read_pfile: genotypes := 'columns' is not compatible with orient := 'genotype' "
				    "(genotype mode already produces scalar output)");
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
			bind_data->tidy_sample_col_to_psam_col.push_back(i);
		}

		bind_data->tidy_sample_col_start = 5;

		// Genotype column (scalar TINYINT or ARRAY(TINYINT,2) when phased)
		names.push_back("genotype");
		return_types.push_back(bind_data->include_phased ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
		                                                 : LogicalType::TINYINT);
		bind_data->tidy_genotype_col = names.size() - 1;
		bind_data->tidy_total_cols = names.size();
	} else if (bind_data->orient_mode == OrientMode::SAMPLE) {
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

		if (bind_data->genotype_mode == GenotypeMode::COLUMNS) {
			// Columns mode: one scalar TINYINT column per effective variant
			if (effective_variant_ct > MAX_GENOTYPE_COLUMNS && !bind_data->has_variant_filter &&
			    !bind_data->region.active) {
				throw InvalidInputException(
				    "read_pfile: genotypes := 'columns' with orient := 'sample' would create %u columns (limit: %u). "
				    "Use variants := [...] or region := '...' to select a subset of variants.",
				    effective_variant_ct, MAX_GENOTYPE_COLUMNS);
			}

			bind_data->columns_mode_first_geno_col = names.size();
			bind_data->columns_mode_geno_col_count = effective_variant_ct;

			// Check for duplicate variant IDs (variant IDs can duplicate unlike sample IIDs)
			std::unordered_set<string> seen_names;
			for (uint32_t ev = 0; ev < effective_variant_ct; ev++) {
				uint32_t vidx = bind_data->has_effective_variant_list ? bind_data->effective_variant_indices[ev] : ev;
				auto id = bind_data->variants.GetId(vidx);
				string col_name;
				if (id.empty()) {
					col_name = bind_data->variants.GetChrom(vidx) + ":" + std::to_string(bind_data->variants.GetPos(vidx));
				} else {
					col_name = id;
				}
				if (!seen_names.insert(col_name).second) {
					throw InvalidInputException(
					    "read_pfile: genotypes := 'columns' with orient := 'sample' requires unique variant "
					    "identifiers, but '%s' appears more than once. Use variants := [...] to select unique variants.",
					    col_name);
				}
				bind_data->genotype_column_names.push_back(col_name);
				names.push_back(col_name);
				return_types.push_back(LogicalType::TINYINT);
			}

			bind_data->sample_orient_total_cols = names.size();
		} else {
			LogicalType sample_elem_type = bind_data->include_phased
			                                   ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
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
		uint64_t matrix_size =
		    static_cast<uint64_t>(effective_variant_ct) * static_cast<uint64_t>(output_sample_ct);
		if (matrix_size > 128ULL * 1024 * 1024) {
			throw InvalidInputException("read_pfile: orient := 'sample' would require %llu genotype values "
			                            "(%u variants x %u samples). Use variants := [...] or samples := [...] to reduce.",
			                            static_cast<unsigned long long>(matrix_size), effective_variant_ct,
			                            output_sample_ct);
		}

		// Init temporary PgenReader for pre-reading
		plink2::PgenFileInfo tmp_pgfi;
		plink2::PreinitPgfi(&tmp_pgfi);

		char errstr_buf2[plink2::kPglErrstrBufBlen];
		plink2::PgenHeaderCtrl header_ctrl2;
		uintptr_t pgfi_alloc_ct2 = 0;

		err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr, bind_data->raw_variant_ct,
		                             bind_data->raw_sample_ct, &header_ctrl2, &tmp_pgfi, &pgfi_alloc_ct2,
		                             errstr_buf2);
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

		// Pre-read genotypes for each effective variant
		bind_data->genotype_matrix.resize(effective_variant_ct);
		for (uint32_t ev = 0; ev < effective_variant_ct; ev++) {
			uint32_t vidx = bind_data->has_effective_variant_list ? bind_data->effective_variant_indices[ev] : ev;

			const uintptr_t *si_ptr = bind_data->has_sample_subset ? preread_subset.SampleInclude() : nullptr;

			if (bind_data->include_phased) {
				uint32_t phasepresent_ct = 0;
				err = plink2::PgrGetP(si_ptr, pssi2, output_sample_ct, vidx, &tmp_pgr,
				                      genovec_buf2.As<uintptr_t>(), preread_phasepresent.As<uintptr_t>(),
				                      preread_phaseinfo.As<uintptr_t>(), &phasepresent_ct);
				if (err != plink2::kPglRetSuccess) {
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&tmp_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&tmp_pgfi, &ce);
					throw IOException("read_pfile: PgrGetP failed for variant %u during sample-orient pre-read",
					                  vidx);
				}
				plink2::GenoarrToBytesMinus9(genovec_buf2.As<uintptr_t>(), output_sample_ct, tmp_bytes.data());
				UnpackPhasedGenotypes(tmp_bytes.data(), preread_phasepresent.As<uintptr_t>(),
				                      preread_phaseinfo.As<uintptr_t>(), output_sample_ct,
				                      preread_phased_pairs.data());
				bind_data->genotype_matrix[ev].assign(preread_phased_pairs.begin(),
				                                      preread_phased_pairs.begin() + output_sample_ct * 2);
			} else {
				err = plink2::PgrGet(si_ptr, pssi2, output_sample_ct, vidx, &tmp_pgr,
				                     genovec_buf2.As<uintptr_t>());
				if (err != plink2::kPglRetSuccess) {
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&tmp_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&tmp_pgfi, &ce);
					throw IOException("read_pfile: PgrGet failed for variant %u during sample-orient pre-read",
					                  vidx);
				}
				plink2::GenoarrToBytesMinus9(genovec_buf2.As<uintptr_t>(), output_sample_ct, tmp_bytes.data());
				bind_data->genotype_matrix[ev].assign(tmp_bytes.begin(),
				                                      tmp_bytes.begin() + output_sample_ct);
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

		if (bind_data->genotype_mode == GenotypeMode::COLUMNS) {
			// Columns mode: one scalar TINYINT column per output sample
			if (output_sample_ct > MAX_GENOTYPE_COLUMNS && !bind_data->has_sample_subset) {
				throw InvalidInputException(
				    "read_pfile: genotypes := 'columns' would create %u columns (limit: %u). "
				    "Use samples := [...] to select a subset of samples.",
				    output_sample_ct, MAX_GENOTYPE_COLUMNS);
			}

			names = {"CHROM", "POS", "ID", "REF", "ALT"};
			return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
			                LogicalType::VARCHAR, LogicalType::VARCHAR};

			bind_data->columns_mode_first_geno_col = names.size();
			bind_data->columns_mode_geno_col_count = output_sample_ct;

			for (uint32_t s = 0; s < output_sample_ct; s++) {
				uint32_t file_idx = bind_data->has_sample_subset ? bind_data->sample_indices[s] : s;
				string col_name = bind_data->sample_info.iids[file_idx];
				bind_data->genotype_column_names.push_back(col_name);
				names.push_back(col_name);
				return_types.push_back(LogicalType::TINYINT);
			}
		} else {
			LogicalType variant_elem_type = bind_data->include_phased
			                                    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
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

	if (bind_data.orient_mode == OrientMode::GENOTYPE) {
		state->total_count = bind_data.EffectiveVariantCt();
		// Check if genotype column is projected
		state->need_genotypes = false;
		for (auto col_id : input.column_ids) {
			if (col_id == bind_data.tidy_genotype_col) {
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
				// In columns mode, any column in the genotype range triggers genotype reading
				if (col_id >= bind_data.columns_mode_first_geno_col &&
				    col_id < bind_data.columns_mode_first_geno_col + bind_data.columns_mode_geno_col_count) {
					state->need_genotypes = true;
					break;
				}
			} else if (col_id == PfileBindData::GENOTYPES_COL) {
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

	if (!gstate.need_genotypes || bind_data.orient_mode == OrientMode::SAMPLE) {
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

			// Read genotype data if needed
			bool genotypes_read = false;
			if (gstate.need_genotypes && lstate.initialized) {
				const uintptr_t *sample_include =
				    bind_data.has_sample_subset ? lstate.sample_include_buf.As<uintptr_t>() : nullptr;

				if (bind_data.include_phased) {
					plink2::PglErr err = plink2::PgrGetP(
					    sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
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
					plink2::PglErr err =
					    plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
					                   lstate.genovec_buf.As<uintptr_t>());
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
							int8_t geno = lstate.genotype_bytes[sample_pos];
							if (geno == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
							}
						} else {
							FlatVector::SetNull(vec, rows_emitted, true);
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
								if (a1 == -9) {
									pair_validity.SetInvalid(pair_idx);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = lstate.phased_pairs[s * 2 + 1];
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
								if (a1 == -9) {
									pair_validity.SetInvalid(pair_idx);
									allele_data[allele_base] = 0;
									allele_data[allele_base + 1] = 0;
								} else {
									allele_data[allele_base] = a1;
									allele_data[allele_base + 1] = lstate.phased_pairs[s * 2 + 1];
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
								if (geno == -9) {
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
								if (geno == -9) {
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
							int8_t geno = lstate.genotype_bytes[sample_pos];
							if (geno == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
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

	if (lstate.tidy_done) {
		output.SetCardinality(0);
		return;
	}

	auto &column_ids = gstate.column_ids;
	uint32_t total_effective_variants = gstate.total_count;
	uint32_t output_sample_ct = bind_data.OutputSampleCt();

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

			if (bind_data.include_phased) {
				plink2::PglErr err = plink2::PgrGetP(
				    sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
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
				plink2::PglErr err =
				    plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
				                   lstate.genovec_buf.As<uintptr_t>());
				if (err != plink2::kPglRetSuccess) {
					throw IOException("read_pfile: PgrGet failed for variant %u", vidx);
				}
				plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
				                             lstate.genotype_bytes.data());
			}
			lstate.tidy_variant_loaded = true;
		}

		// Emit rows for samples within current variant
		while (lstate.tidy_current_sample < output_sample_ct && rows_emitted < STANDARD_VECTOR_SIZE) {
			// Map scan-order sample index to file-order sample index.
			// When subsetting, sample_indices is sorted to match pgenlib's
			// ascending-bit-order output, so genotype_bytes[i] corresponds
			// to sample_indices[i].
			uint32_t sample_file_idx = bind_data.has_sample_subset
			                               ? bind_data.sample_indices[lstate.tidy_current_sample]
			                               : lstate.tidy_current_sample;

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
				} else if (file_col == bind_data.tidy_genotype_col) {
					// Genotype column
					if (gstate.need_genotypes && lstate.tidy_variant_loaded) {
						if (bind_data.include_phased) {
							// ARRAY(TINYINT, 2) per row
							auto &allele_vec = ArrayVector::GetEntry(vec);
							auto *allele_data = FlatVector::GetData<int8_t>(allele_vec);
							idx_t allele_base = rows_emitted * 2;
							int8_t a1 = lstate.phased_pairs[lstate.tidy_current_sample * 2];
							if (a1 == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
								allele_data[allele_base] = 0;
								allele_data[allele_base + 1] = 0;
							} else {
								allele_data[allele_base] = a1;
								allele_data[allele_base + 1] =
								    lstate.phased_pairs[lstate.tidy_current_sample * 2 + 1];
							}
						} else {
							// Scalar TINYINT
							int8_t geno = lstate.genotype_bytes[lstate.tidy_current_sample];
							if (geno == -9) {
								FlatVector::SetNull(vec, rows_emitted, true);
							} else {
								FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
							}
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
			uint32_t sample_file_idx =
			    bind_data.has_sample_subset ? bind_data.sample_indices[sample_pos] : sample_pos;

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
						int8_t geno = bind_data.genotype_matrix[variant_pos][sample_pos];
						if (geno == -9) {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<int8_t>(vec)[rows_emitted] = geno;
						}
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
				} else if (bind_data.genotype_mode != GenotypeMode::COLUMNS &&
				           file_col == bind_data.sample_orient_genotypes_col) {
					// Genotypes column — ARRAY or LIST of genotypes across variants
					if (!gstate.need_genotypes) {
						FlatVector::SetNull(vec, rows_emitted, true);
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
									allele_data[allele_base + 1] =
									    bind_data.genotype_matrix[v][sample_pos * 2 + 1];
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
									allele_data[allele_base + 1] =
									    bind_data.genotype_matrix[v][sample_pos * 2 + 1];
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
	read_pfile.named_parameters["tidy"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["orient"] = LogicalType::VARCHAR;
	read_pfile.named_parameters["dosages"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["phased"] = LogicalType::BOOLEAN;
	read_pfile.named_parameters["region"] = LogicalType::VARCHAR;
	// Accept ANY for samples and variants — type dispatch handled in bind
	read_pfile.named_parameters["samples"] = LogicalType::ANY;
	read_pfile.named_parameters["variants"] = LogicalType::ANY;
	read_pfile.named_parameters["genotypes"] = LogicalType::VARCHAR;

	loader.RegisterFunction(read_pfile);
}

} // namespace duckdb
