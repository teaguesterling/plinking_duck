#include "pfile_reader.hpp"
#include "duckdb_compat.hpp"
#include "plink_common.hpp"
#include "plink_profile.hpp"
#include "pvar_reader.hpp"
#include "psam_reader.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/function_set.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/materialized_query_result.hpp"

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
	//! All data rows from the .psam/.fam, in file order (TEXT path only).
	//! Each inner vector has one field per column. Empty for the parquet path, whose
	//! column values are loaded projection-aware at InitGlobal (PfileGlobalState).
	vector<vector<string>> rows;

	//! Source path for the parquet/non-native path — the projected column data is
	//! (re-)queried from here at InitGlobal, not materialized at bind. Empty for text.
	string source_path;

	//! IID / FID column positions in the header (for the parquet path).
	idx_t iid_col_idx = DConstants::INVALID_INDEX;
	idx_t fid_col_idx = DConstants::INVALID_INDEX;

	//! Index of SEX column in the header (DConstants::INVALID_INDEX if absent)
	idx_t sex_col_idx = DConstants::INVALID_INDEX;
	//! Indices of PAT/MAT columns for special missing-value handling
	vector<idx_t> parent_col_indices;

	//! True when the psam came from a parquet/non-native source: column data is not
	//! in `rows`; the scan reads it from the projection-aware CDC in the global state.
	bool from_parquet = false;
};

//! Missing value sentinels (same as psam_reader)
static bool PfileIsMissingValue(const string &val) {
	return val.empty() || val == "." || val == "NA" || val == "na";
}

//! Load sample metadata SCHEMA (not data) from a non-native source (parquet/csv/table)
//! for genotype/sample orient modes. The column values are loaded projection-aware at
//! InitGlobal (see BuildProjectedPsamCdc), so a wide psam never materializes columns
//! the query does not select. Only the header, special-column indices, and sample
//! count are established here.
static PfileSampleMetadata LoadPfileSampleMetadataFromSource(ClientContext &context, const string &source,
                                                             SampleInfo &sample_info_out) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto escaped = StringUtil::Replace(source, "'", "''");

	// Schema only — LIMIT 0 reads the footer/metadata, no row groups.
	auto schema_res = conn.Query("SELECT * FROM '" + escaped + "' LIMIT 0");
	if (schema_res->HasError()) {
		throw IOException("read_pfile: failed to read source '%s': %s", source, schema_res->GetError());
	}

	PfileSampleMetadata meta;
	meta.from_parquet = true;
	meta.source_path = source;
	meta.header.format = PsamFormat::PSAM_IID;
	meta.header.column_names = schema_res->names;
	meta.header.column_types = schema_res->types;

	idx_t iid_idx = DConstants::INVALID_INDEX;
	for (idx_t i = 0; i < schema_res->names.size(); i++) {
		auto lower = StringUtil::Lower(schema_res->names[i]);
		if (lower == "sex") {
			meta.sex_col_idx = i;
		} else if (lower == "pat" || lower == "mat") {
			meta.parent_col_indices.push_back(i);
		}
		if (lower == "iid") {
			iid_idx = i;
		} else if (lower == "fid") {
			meta.fid_col_idx = i;
		}
	}
	if (iid_idx == DConstants::INVALID_INDEX) {
		throw IOException("read_pfile: source '%s' has no IID column (found: %s)", source,
		                  StringUtil::Join(schema_res->names, ", "));
	}
	meta.iid_col_idx = iid_idx;

	auto cnt_res = conn.Query("SELECT COUNT(*) FROM '" + escaped + "'");
	if (cnt_res->HasError()) {
		throw IOException("read_pfile: failed to count rows in source '%s': %s", source, cnt_res->GetError());
	}
	sample_info_out.sample_ct = cnt_res->GetValue(0, 0).GetValue<int64_t>();
	// iids are NOT materialized here — see EnsureSourceIids (only a VARCHAR sample
	// subset needs them; IID output reads from the projected CDC at scan).
	return meta;
}

//! Extract the string_t values of one VARCHAR-castable column of a fetched/scanned
//! chunk into a std::string vector (NULL/invalid -> empty string).
static void ExtractStringColumn(ClientContext &context, DataChunk &chunk, idx_t chunk_col, idx_t n,
                                vector<string> &out) {
	Vector varchar_vec(LogicalType::VARCHAR);
	Vector *src = &chunk.data[chunk_col];
	if (src->GetType().id() != LogicalTypeId::VARCHAR) {
		VectorOperations::Cast(context, *src, varchar_vec, n);
		src = &varchar_vec;
	}
	UnifiedVectorFormat fmt;
	src->ToUnifiedFormat(n, fmt);
	auto data = UnifiedVectorFormat::GetData<string_t>(fmt);
	for (idx_t row = 0; row < n; row++) {
		idx_t vidx = fmt.sel->get_index(row);
		out.push_back(fmt.validity.RowIsValid(vidx) ? data[vidx].GetString() : string());
	}
}

//! Lazily materialize SampleInfo::iids (and fids) for a parquet-source psam, by
//! querying just the IID (and FID) columns. No-op if already populated or not a
//! parquet source. Only a `samples := [IID strings]` subset needs this — IID/FID
//! *output* reads natively from the projected CDC and never triggers it.
static void EnsureSourceIids(ClientContext &context, const PfileSampleMetadata &meta, SampleInfo &sample_info) {
	if (!sample_info.iids.empty() || !meta.from_parquet) {
		return;
	}
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto escaped = StringUtil::Replace(meta.source_path, "'", "''");
	const bool has_fid = meta.fid_col_idx != DConstants::INVALID_INDEX;
	auto quote = [](const string &n) {
		return "\"" + StringUtil::Replace(n, "\"", "\"\"") + "\"";
	};
	string cols = quote(meta.header.column_names[meta.iid_col_idx]);
	if (has_fid) {
		cols += ", " + quote(meta.header.column_names[meta.fid_col_idx]);
	}
	auto result = conn.Query("SELECT " + cols + " FROM '" + escaped + "'");
	if (result->HasError()) {
		throw IOException("read_pfile: failed to load IIDs from source '%s': %s", meta.source_path, result->GetError());
	}
	sample_info.iids.reserve(sample_info.sample_ct);
	if (has_fid) {
		sample_info.fids.reserve(sample_info.sample_ct);
	}
	unique_ptr<DataChunk> chunk;
	while ((chunk = result->Fetch()) != nullptr && chunk->size() > 0) {
		const idx_t n = chunk->size();
		ExtractStringColumn(context, *chunk, 0, n, sample_info.iids);
		if (has_fid) {
			ExtractStringColumn(context, *chunk, 1, n, sample_info.fids);
		}
	}
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

//! One variant-sharded pfile source. A single-file read has exactly one; a
//! LIST(VARCHAR) read has one per prefix, row-concatenated in list order. All
//! sources share an identical .psam by contract (same samples, same order) — only
//! the variants differ. Sample metadata lives on PfileBindData (from sources[0]).
struct PfileSource {
	string pgen_path;
	string pvar_path;
	string psam_path;

	// This shard's variant metadata (post region/variant filter is applied via the
	// effective list below; `variants` itself holds the full or region-loaded index).
	VariantMetadataIndex variants;

	// pgenlib header: this file's variant count. raw_sample_ct is shared (PfileBindData).
	uint32_t raw_variant_ct = 0;

	// Effective variant list for this source (intersection of region + variant filter).
	// When has_effective_variant_list is false, scan all raw_variant_ct sequentially.
	bool has_effective_variant_list = false;
	vector<uint32_t> effective_variant_indices;

	//! Number of effective variants in this source (after filtering).
	uint32_t EffectiveVariantCt() const {
		return has_effective_variant_list ? static_cast<uint32_t>(effective_variant_indices.size()) : raw_variant_ct;
	}

	//! Map a 0-based effective position within this source to its pgen variant index.
	uint32_t ResolveVariantIdx(uint32_t effective_pos) const {
		return has_effective_variant_list ? effective_variant_indices[effective_pos] : effective_pos;
	}
};

// How the sample sets of multiple sources are combined. Only IMPLICIT and
// IDENTICAL are implemented; the others are reserved for real sample-set joins.
//   IMPLICIT    - trust the caller: same samples, same order, positional (default,
//                 matches the existing multi-file contract; no per-shard psam I/O).
//   IDENTICAL   - verify every shard's IIDs are equal (same values, same order),
//                 then combine positionally. A mismatch is a hard error.
//   UNION       - outer-join samples across shards (not yet implemented).
//   INTERSECT   - inner-join samples across shards (not yet implemented).
//   CONCATENATE - stack DIFFERENT samples (sample-sharded) (not yet implemented).
enum class CombineSamplesMode : uint8_t { IMPLICIT, IDENTICAL, UNION, INTERSECT, CONCATENATE };

static CombineSamplesMode ResolveCombineSamplesMode(const string &s) {
	string v = StringUtil::Lower(s);
	StringUtil::Trim(v);
	if (v == "implicit") {
		return CombineSamplesMode::IMPLICIT;
	}
	if (v == "identical") {
		return CombineSamplesMode::IDENTICAL;
	}
	if (v == "union") {
		return CombineSamplesMode::UNION;
	}
	if (v == "intersect") {
		return CombineSamplesMode::INTERSECT;
	}
	if (v == "concatenate") {
		return CombineSamplesMode::CONCATENATE;
	}
	throw InvalidInputException("read_pfile: unknown combine_samples := '%s' (expected 'implicit', 'identical', "
	                            "'union', 'intersect', or 'concatenate')",
	                            s);
}

struct PfileBindData : public TableFunctionData {
	// One or more variant-sharded sources, row-concatenated in list order.
	// Single-file reads have exactly one; sources[0] is authoritative for the
	// shared sample metadata (identical .psam across shards by contract).
	vector<PfileSource> sources;

	// Cumulative effective-variant offsets across sources (size sources.size()+1).
	// Global effective position g in [variant_offsets[i], variant_offsets[i+1]) maps
	// to source i, local position g - variant_offsets[i]. Used by multi-file scan.
	vector<uint32_t> variant_offsets;

	// Sample metadata (basic: for IID lookups and default mode) — from sources[0]
	SampleInfo sample_info;

	// Full sample metadata (for genotype orient mode column output) — from sources[0]
	PfileSampleMetadata sample_metadata;

	// pgenlib header: sample count (shared — identical across all sources by contract)
	uint32_t raw_sample_ct = 0;

	// Mode
	OrientMode orient_mode = OrientMode::VARIANT;

	// How multiple sources' sample sets combine (multi-file only).
	CombineSamplesMode combine_samples = CombineSamplesMode::IMPLICIT;

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

	// Variant filtering — the `variants` named param, resolved per source in bind.
	// The effective variant list (region ∩ variants) lives on each PfileSource.
	bool has_variant_filter = false;

	// Count-based filtering (af_range, ac_range)
	CountFilter count_filter;
	unique_ptr<SampleSubset> count_filter_subset; // shared bind-time subset for PgrGetCounts

	// Genotype range filtering (genotype_range)
	GenotypeRangeFilter genotype_filter;

	// Sample-orient mode: pre-read genotype matrix (variant × sample)
	// genotype_matrix[effective_vidx][sample_idx] = genotype value (-9 = missing)
	vector<vector<int8_t>> genotype_matrix;
	// Dosage variant of genotype_matrix (-9.0 = missing)
	vector<vector<double>> dosage_matrix;

	// Sample-orient row-level genotype_range filter: output sample positions that
	// survive the filter (a sample is kept iff at least one of its genotypes
	// satisfies genotype_range, or is missing when include_missing is set).
	// When has_sample_keep is false, all OutputSampleCt() samples are emitted.
	bool has_sample_keep = false;
	vector<uint32_t> sample_keep;

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

	//! Total effective variants across all sources (global scan count). For a
	//! single-file read this equals sources[0].EffectiveVariantCt().
	uint32_t EffectiveVariantCt() const {
		return variant_offsets.empty() ? 0 : variant_offsets.back();
	}

	//! Convenience accessor for the primary (and, single-file, only) source.
	PfileSource &Primary() {
		return sources[0];
	}
	const PfileSource &Primary() const {
		return sources[0];
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

//! A unit of scan work bounded to a single source (multi-file variant/genotype
//! orient). local_start/local_end are 0-based effective positions within that
//! source. Threads claim batches atomically; a batch never spans a file, so a
//! thread reopens the pgen reader at most once per claimed batch.
struct ScanBatch {
	uint32_t source_idx;
	uint32_t local_start;
	uint32_t local_end;
};

struct PfileGlobalState : public GlobalTableFunctionState {
	// For variant/sample-orient mode: atomic counter
	std::atomic<uint32_t> next_idx {0};
	uint32_t total_count = 0; // variants (variant/genotype mode) or samples (sample mode)

	// Multi-file (sources.size() > 1) variant/genotype orient: precomputed batches,
	// each bounded to one source. Claimed via next_idx (as a batch index). Empty for
	// single-file reads, which keep their original claim loops.
	vector<ScanBatch> batches;

	// Projection
	bool need_genotypes = false;
	bool need_pgen_reader = false; // genotypes OR count filter
	vector<column_t> column_ids;
	OrientMode orient_mode = OrientMode::VARIANT;
	uint32_t max_threads_config = 0;

	// Projection-aware psam column data (parquet source only). Built once here with
	// ONLY the psam columns the query projects, so a wide biobank psam never
	// materializes/reads unused covariate columns. Read at scan via a per-thread
	// cached FetchChunk (thread-safe: FetchChunk on a shared const CDC).
	unique_ptr<ColumnDataCollection> psam_cdc;
	vector<idx_t> psam_chunk_start_row; // cumulative chunk offsets (size ChunkCount()+1)
	vector<idx_t> psam_col_to_cdc;      // header col idx -> psam_cdc col idx (INVALID if not projected)

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

	// Multi-file: index of the source this thread's pgr/pgfi are currently open on
	// (DConstants::INVALID_INDEX = none open yet). Reopened on a source boundary.
	idx_t current_source_idx = DConstants::INVALID_INDEX;

	// Multi-file batch-claim position (variant + genotype orient). mf_batch is the
	// claimed index into PfileGlobalState::batches; mf_local is the next effective
	// position within that batch's [local_start, local_end). mf_need_claim triggers
	// claiming a fresh batch.
	uint32_t mf_batch = 0;
	uint32_t mf_local = 0;
	bool mf_need_claim = true;

	// Genotype-orient batch claiming state
	uint32_t batch_start = 0;               // first effective variant index in current batch
	uint32_t batch_end = 0;                 // one-past-end effective variant index in current batch
	uint32_t current_variant_in_batch = 0;  // current effective variant index within batch
	uint32_t current_sample_in_variant = 0; // sample index within current variant
	bool batch_variant_loaded = false;      // genotypes decoded for current variant in batch
	bool geno_range_all_pass = true;        // per-variant flag for genotype_range optimization
	bool batch_exhausted = true;            // true initially to trigger first batch claim

	// Per-thread cache of one chunk of the projected psam CDC (PfileGlobalState::psam_cdc),
	// for reading psam column values at scan time. Populated lazily via FetchChunk.
	DataChunk psam_chunk;
	idx_t psam_chunk_idx = DConstants::INVALID_INDEX;

	~PfileLocalState() {
		if (initialized) {
			plink2::PglErr reterr = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr, &reterr);
			plink2::CleanupPgfi(&pgfi, &reterr);
		}
	}
};

// ---------------------------------------------------------------------------
// Per-source loading (companion discovery + pgenlib header + variant metadata)
// ---------------------------------------------------------------------------

//! Discover companion files, read the .pgen header, and load variant metadata for
//! ONE source (shard). `override_*` supply explicit paths (single-file only; empty
//! for shards in a list). `need_psam` gates .psam discovery — only the first source
//! needs it (shared sample metadata). Outputs this file's sample count via
//! `raw_sample_ct_out` for the caller's cross-source equality check.
static PfileSource LoadPfileSource(ClientContext &context, FileSystem &fs, const string &prefix,
                                   const string &override_pgen, const string &override_pvar,
                                   const string &override_psam, bool need_psam, const RegionFilter &region,
                                   uint32_t &raw_sample_ct_out) {
	PfileSource src;
	// Resolve explicit overrides against file_search_path too (keep the literal
	// if not found, so downstream open produces the natural error message).
	src.pgen_path = override_pgen;
	src.pvar_path = override_pvar;
	src.psam_path = override_psam;
	for (auto *p : {&src.pgen_path, &src.pvar_path, &src.psam_path}) {
		if (!p->empty()) {
			auto resolved = ResolveExistingPath(context, fs, *p);
			if (!resolved.empty()) {
				*p = resolved;
			}
		}
	}

	// Discover .pgen from prefix if not explicitly provided. Resolve the base
	// against file_search_path ONCE (mirrors read_csv), then derive all
	// companions from the resolved concrete base so they share its directory —
	// never search-resolve each companion independently (that could mix files
	// from different search dirs).
	string eff_prefix = prefix;
	if (src.pgen_path.empty() && !prefix.empty()) {
		string resolved = ResolveExistingPath(context, fs, prefix + ".pgen");
		if (!resolved.empty()) {
			src.pgen_path = resolved;
			eff_prefix = resolved.substr(0, resolved.size() - 5); // strip ".pgen"
		} else {
			resolved = ResolveExistingPath(context, fs, prefix);
			if (!resolved.empty()) {
				src.pgen_path = resolved;
				eff_prefix = resolved;
			} else {
				throw InvalidInputException("read_pfile: cannot find .pgen file for prefix '%s' (tried '%s')", prefix,
				                            prefix + ".pgen");
			}
		}
	}
	if (src.pgen_path.empty()) {
		throw InvalidInputException("read_pfile: no .pgen file path provided");
	}

	// Discover .pvar/.bim
	if (src.pvar_path.empty()) {
		if (!eff_prefix.empty()) {
			src.pvar_path = FindCompanionFileWithParquet(context, fs, eff_prefix + ".pgen", {".pvar", ".bim"});
			if (src.pvar_path.empty()) {
				for (auto &ext : {".pvar", ".bim"}) {
					auto candidate = eff_prefix + ext;
					if (fs.FileExists(candidate)) {
						src.pvar_path = candidate;
						break;
					}
				}
			}
		}
		if (src.pvar_path.empty()) {
			src.pvar_path = FindCompanionFileWithParquet(context, fs, src.pgen_path, {".pvar", ".bim"});
		}
		if (src.pvar_path.empty()) {
			throw InvalidInputException("read_pfile: cannot find .pvar or .bim file for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            prefix.empty() ? src.pgen_path : prefix);
		}
	}

	// Discover .psam/.fam (only the first source needs it — shared sample metadata)
	if (need_psam && src.psam_path.empty()) {
		if (!eff_prefix.empty()) {
			src.psam_path = FindCompanionFileWithParquet(context, fs, eff_prefix + ".pgen", {".psam", ".fam"});
		}
		// Fall back to deriving from the resolved .pgen path. This is required
		// when eff_prefix is itself a full .pgen path (e.g. from a '*.pgen' glob):
		// eff_prefix + ".pgen" would double the extension, so ReplaceExtension on
		// src.pgen_path is what finds the companion. Mirrors the .pvar block.
		if (src.psam_path.empty()) {
			src.psam_path = FindCompanionFileWithParquet(context, fs, src.pgen_path, {".psam", ".fam"});
		}
		if (src.psam_path.empty()) {
			throw InvalidInputException("read_pfile: cannot find .psam or .fam file for '%s' "
			                            "(use psam := 'path' to specify explicitly)",
			                            prefix.empty() ? src.pgen_path : prefix);
		}
	}

	// --- Read the .pgen header (counts) ---
	plink2::PgenFileInfo pgfi;
	plink2::PreinitPgfi(&pgfi);
	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err = plink2::PgfiInitPhase1(src.pgen_path.c_str(), nullptr, UINT32_MAX, UINT32_MAX, &header_ctrl,
	                                            &pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("read_pfile: failed to open '%s': %s", src.pgen_path, errstr_buf);
	}

	src.raw_variant_ct = pgfi.raw_variant_ct;
	raw_sample_ct_out = pgfi.raw_sample_ct;

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
		throw IOException("read_pfile: failed to initialize '%s' (phase 2): %s", src.pgen_path, errstr_buf);
	}

	// --- Load variant metadata (region pushdown for parquet) ---
	if (region.active && IsParquetFile(src.pvar_path)) {
		src.variants =
		    LoadVariantMetadataFromParquetRegion(context, src.pvar_path, region.chrom, region.start, region.end,
		                                         static_cast<idx_t>(src.raw_variant_ct), "read_pfile");
	} else {
		src.variants = LoadVariantMetadata(context, src.pvar_path, "read_pfile");
	}

	if (src.variants.variant_ct != src.raw_variant_ct) {
		throw InvalidInputException("read_pfile: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            src.raw_variant_ct, src.pvar_path,
		                            static_cast<unsigned long long>(src.variants.variant_ct));
	}

	return src;
}

//! Build one source's effective variant list (intersection of region + variants
//! filter). `variant_set` (when `use_variant_set`) holds the resolved `variants`
//! param indices FOR THIS SOURCE. Mirrors the single-file logic, keyed on `src`.
static void BuildEffectiveVariantList(PfileSource &src, const RegionFilter &region,
                                      const std::unordered_set<uint32_t> &variant_set, bool use_variant_set) {
	bool any_filter_active = region.active || use_variant_set;
	if (!any_filter_active) {
		return; // scan all raw_variant_ct sequentially
	}
	src.has_effective_variant_list = true;

	// Case A: sparse/subset pvar index (parquet region pushdown). Use IsDense()
	// (not vidx_map.empty()) so a zero-match region — an empty subset — takes this
	// path and yields an empty effective list, instead of falling through to the
	// dense Case C, which would index the empty metadata vectors out of bounds.
	if (!src.variants.IsDense()) {
		src.effective_variant_indices.reserve(src.variants.vidx_map.size());
		for (auto &kv : src.variants.vidx_map) {
			if (use_variant_set && variant_set.find(kv.first) == variant_set.end()) {
				continue;
			}
			src.effective_variant_indices.push_back(kv.first);
		}
		std::sort(src.effective_variant_indices.begin(), src.effective_variant_indices.end());
	} else if (region.active && !src.variants.chrom_offsets.empty()) {
		// Case B: dense + region + chrom_offsets → O(log N) binary-search bounds.
		auto it = src.variants.chrom_offsets.find(region.chrom);
		if (it != src.variants.chrom_offsets.end()) {
			idx_t lo_local = it->second.first;
			idx_t hi_local = it->second.second;
			auto &positions = src.variants.positions;
			idx_t lo = lo_local, hi = hi_local;
			while (lo < hi) {
				idx_t mid = lo + (hi - lo) / 2;
				if (positions[mid] < region.start) {
					lo = mid + 1;
				} else {
					hi = mid;
				}
			}
			idx_t start_idx = lo;
			lo = lo_local;
			hi = hi_local;
			while (lo < hi) {
				idx_t mid = lo + (hi - lo) / 2;
				if (positions[mid] <= region.end) {
					lo = mid + 1;
				} else {
					hi = mid;
				}
			}
			idx_t end_idx = lo;
			src.effective_variant_indices.reserve(end_idx - start_idx);
			for (idx_t vidx = start_idx; vidx < end_idx; vidx++) {
				if (use_variant_set && variant_set.find(static_cast<uint32_t>(vidx)) == variant_set.end()) {
					continue;
				}
				src.effective_variant_indices.push_back(static_cast<uint32_t>(vidx));
			}
		}
	} else {
		// Case C: variant filter only (no region), or region without chrom_offsets.
		for (uint32_t vidx = 0; vidx < src.raw_variant_ct; vidx++) {
			if (region.active) {
				if (src.variants.GetChrom(vidx) != region.chrom) {
					continue;
				}
				int64_t pos = src.variants.GetPos(vidx);
				if (pos < region.start || pos > region.end) {
					continue;
				}
			}
			if (use_variant_set && variant_set.find(vidx) == variant_set.end()) {
				continue;
			}
			src.effective_variant_indices.push_back(vidx);
		}
	}
}

// ---------------------------------------------------------------------------
// Bind function
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PfileBind(ClientContext &context, TableFunctionBindInput &input,
                                          vector<LogicalType> &return_types, vector<string> &names) {
	BindPhaseTimer bind_timer("PfileBind(total)");
	auto bind_data = make_uniq<PfileBindData>();

	// --- Resolve file path(s) ---
	// First positional argument is a prefix (VARCHAR) or a list of prefixes
	// (LIST(VARCHAR), row-concatenated in order). Empty inputs → only named params.
	auto &fs = FileSystem::GetFileSystem(context);

	vector<string> prefixes;
	if (!input.inputs.empty()) {
		prefixes = ResolvePathList(input.inputs[0], "read_pfile");
		// Expand glob patterns and protocol URLs (e.g. pathmacro:) to concrete
		// local prefixes/paths, so a single URL/glob can fan out to a sharded
		// fileset (mirrors read_csv). Plain prefixes pass through unchanged.
		prefixes = ExpandPathInputs(context, fs, prefixes, "read_pfile");
	} else {
		prefixes.push_back("");
	}
	const bool multi_file = prefixes.size() > 1;

	// Named parameters can override individual paths (single-file only)
	string orient_str;
	string override_pgen, override_pvar, override_psam;
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pgen") {
			override_pgen = kv.second.GetValue<string>();
		} else if (kv.first == "pvar") {
			override_pvar = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			override_psam = kv.second.GetValue<string>();
		} else if (kv.first == "orient") {
			orient_str = kv.second.GetValue<string>();
		} else if (kv.first == "dosages") {
			bind_data->include_dosages = kv.second.GetValue<bool>();
		} else if (kv.first == "phased") {
			bind_data->include_phased = kv.second.GetValue<bool>();
		} else if (kv.first == "region") {
			bind_data->region = ParseRegion(kv.second.GetValue<string>());
		} else if (kv.first == "combine_samples") {
			bind_data->combine_samples = ResolveCombineSamplesMode(kv.second.GetValue<string>());
		}
		// samples, variants, genotypes, af_range, ac_range handled after pgenlib init
	}

	bind_data->orient_mode = ResolveOrientMode(orient_str, "read_pfile");

	// Only IMPLICIT and IDENTICAL are implemented; the join modes are reserved.
	if (bind_data->combine_samples == CombineSamplesMode::UNION ||
	    bind_data->combine_samples == CombineSamplesMode::INTERSECT ||
	    bind_data->combine_samples == CombineSamplesMode::CONCATENATE) {
		throw InvalidInputException("read_pfile: combine_samples := '%s' is not yet implemented "
		                            "(only 'implicit' and 'identical' are supported)",
		                            StringUtil::Lower(input.named_parameters.at("combine_samples").GetValue<string>()));
	}

	if (bind_data->include_dosages && bind_data->include_phased) {
		throw InvalidInputException("read_pfile: dosages and phased cannot both be true");
	}

	// --- Multi-file guards (throw before loading anything) ---
	if (multi_file) {
		// pgen/pvar overrides can't disambiguate across shards: the list elements
		// ARE the per-shard pgens, and each shard has its own (different) variants.
		// A psam override IS allowed — all shards share one sample set by contract,
		// so a single .psam legitimately applies to every shard ("use this psam for
		// all pgens"); it becomes the shared sample metadata (see below).
		if (!override_pgen.empty() || !override_pvar.empty()) {
			throw InvalidInputException(
			    "read_pfile: pgen/pvar overrides cannot be combined with a multi-file list (pgen is the per-shard "
			    "input and pvar differs per shard). A psam override is allowed — it applies to every shard.");
		}
		// `variants` resolves against a single file's variant index; across a
		// multi-file list neither by-ID (an ID lives in only one shard) nor by-index
		// (indices are per-shard, not global) resolves correctly yet. Guard it rather
		// than throw-on-missing or return per-shard duplicates. Use `region :=` for
		// range selection across files. Global `variants` resolution is a follow-up.
		if (input.named_parameters.find("variants") != input.named_parameters.end()) {
			throw InvalidInputException("read_pfile: variants := [...] with a multi-file list is not yet supported; "
			                            "use region := for selection across files");
		}
	}

	// Resolve the `variants` param against the variant index below. Multi-file +
	// `variants` is rejected above, so here it only ever applies to a single source.
	auto variants_it = input.named_parameters.find("variants");
	bind_data->has_variant_filter = variants_it != input.named_parameters.end();

	// --- Build sources (one per prefix) ---
	bind_timer.Note("resolving %llu source(s)", (unsigned long long)prefixes.size());
	for (idx_t si = 0; si < prefixes.size(); si++) {
		uint32_t src_sample_ct = 0;
		// combine_samples := 'identical' needs each shard's .psam discovered so we
		// can compare IIDs below; otherwise only source 0's psam is needed. A psam
		// override supplies one shared .psam for all shards, so there is nothing to
		// discover or compare per shard.
		bool need_psam =
		    (si == 0) || (bind_data->combine_samples == CombineSamplesMode::IDENTICAL && override_psam.empty());
		PfileSource src = LoadPfileSource(context, fs, prefixes[si], si == 0 ? override_pgen : string(),
		                                  si == 0 ? override_pvar : string(), si == 0 ? override_psam : string(),
		                                  need_psam, bind_data->region, src_sample_ct);

		// Free safety check: every shard's .pgen header reports its sample_ct. All
		// shards share an identical .psam by contract; a genuine mismatch would
		// silently overflow the per-thread sample buffers — reject it clearly.
		if (si == 0) {
			bind_data->raw_sample_ct = src_sample_ct;
		} else if (src_sample_ct != bind_data->raw_sample_ct) {
			throw InvalidInputException(
			    "read_pfile: sample count mismatch across files: '%s' has %u samples, but expected %u "
			    "(from '%s'). All listed pfiles must share the same samples.",
			    src.pgen_path, src_sample_ct, bind_data->raw_sample_ct, bind_data->sources[0].pgen_path);
		}

		// Per-source effective variant list (region ∩ variants param, resolved per shard)
		std::unordered_set<uint32_t> variant_set;
		if (bind_data->has_variant_filter) {
			auto resolved =
			    ResolveVariantsParameter(variants_it->second, src.variants, src.raw_variant_ct, "read_pfile");
			for (auto idx : resolved) {
				variant_set.insert(idx);
			}
		}
		BuildEffectiveVariantList(src, bind_data->region, variant_set, bind_data->has_variant_filter);

		bind_data->sources.push_back(std::move(src));
	}

	// Cumulative effective-variant offsets across sources (global scan positions)
	bind_data->variant_offsets.push_back(0);
	{
		uint32_t cum = 0;
		for (auto &s : bind_data->sources) {
			cum += s.EffectiveVariantCt();
			bind_data->variant_offsets.push_back(cum);
		}
	}
	bind_timer.Note("sources built (raw_sample_ct=%u, total effective variants=%u)", bind_data->raw_sample_ct,
	                bind_data->EffectiveVariantCt());

	// --- Load sample info (from sources[0] only — identical psam by contract) ---
	const string &psam_path = bind_data->Primary().psam_path;
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
		bind_data->sample_metadata = LoadPfileSampleMetadata(context, psam_path, bind_data->sample_info);
	} else if (needs_iids) {
		bind_data->sample_info = LoadSampleMetadata(context, psam_path);
	} else {
		bind_data->sample_info = LoadSampleCount(context, psam_path);
	}
	bind_timer.Note("sample metadata loaded (%llu samples, iids=%llu)",
	                (unsigned long long)bind_data->sample_info.sample_ct,
	                (unsigned long long)bind_data->sample_info.iids.size());

	if (bind_data->sample_info.sample_ct != bind_data->raw_sample_ct) {
		throw InvalidInputException("read_pfile: sample count mismatch: .pgen has %u samples, "
		                            ".psam/.fam '%s' has %llu samples",
		                            bind_data->raw_sample_ct, psam_path,
		                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
	}

	// combine_samples := 'identical': verify every shard carries the SAME IIDs in
	// the SAME order as source 0, then combine positionally. Equal .pgen sample_ct
	// (checked above) does not imply equal IIDs; without this a mismatch would
	// silently fuse different individuals into one row (esp. in sample orient).
	// Cost: one .psam parse per shard — the price of validation; 'implicit' (the
	// default) skips it and trusts the caller.
	if (multi_file && bind_data->combine_samples == CombineSamplesMode::IDENTICAL && override_psam.empty()) {
		vector<string> ref_iids_owned;
		const vector<string> *ref_iids = &bind_data->sample_info.iids;
		if (ref_iids->empty()) { // variant orient loads count-only — fetch IIDs for the check
			ref_iids_owned = LoadSampleMetadata(context, bind_data->sources[0].psam_path).iids;
			ref_iids = &ref_iids_owned;
		}
		for (idx_t si = 1; si < bind_data->sources.size(); si++) {
			auto iids = LoadSampleMetadata(context, bind_data->sources[si].psam_path).iids;
			if (iids.size() != ref_iids->size()) {
				throw InvalidInputException(
				    "read_pfile: combine_samples := 'identical' but '%s' has %llu samples vs %llu in '%s'",
				    bind_data->sources[si].psam_path, static_cast<unsigned long long>(iids.size()),
				    static_cast<unsigned long long>(ref_iids->size()), bind_data->sources[0].psam_path);
			}
			for (idx_t k = 0; k < iids.size(); k++) {
				if (iids[k] != (*ref_iids)[k]) {
					throw InvalidInputException(
					    "read_pfile: combine_samples := 'identical' but sample %llu differs: '%s' in '%s' vs '%s' in "
					    "'%s'. All shards must share the same IIDs in the same order (or use combine_samples := "
					    "'implicit' to trust positional alignment).",
					    static_cast<unsigned long long>(k), iids[k], bind_data->sources[si].psam_path, (*ref_iids)[k],
					    bind_data->sources[0].psam_path);
				}
			}
		}
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
			// Lazily build iid_to_idx only when a VARCHAR sample filter is used. On the
			// parquet path the IID strings are not materialized at load, so pull just the
			// IID/FID columns from the source first (no-op on the text path).
			EnsureSourceIids(context, bind_data->sample_metadata, bind_data->sample_info);
			bind_data->sample_info.EnsureIidMap("read_pfile: " + psam_path);
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

	// (variants param + effective variant list were resolved per-source above.)

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

	// --- Parse genotype filter: include_genotypes (canonical) or genotype_range (alias) ---
	{
		auto ig_it = input.named_parameters.find("include_genotypes");
		auto gr_it = input.named_parameters.find("genotype_range");
		bool has_ig = ig_it != input.named_parameters.end();
		bool has_gr = gr_it != input.named_parameters.end();
		if (has_ig && has_gr) {
			throw InvalidInputException(
			    "read_pfile: specify only one of include_genotypes or genotype_range (genotype_range is the numeric "
			    "alias of include_genotypes)");
		}
		if (has_ig || has_gr) {
			if (bind_data->include_dosages) {
				throw InvalidInputException("read_pfile: %s is incompatible with dosages := true",
				                            has_ig ? "include_genotypes" : "genotype_range");
			}
		}
		if (has_ig) {
			ParseIncludeGenotypes(ig_it->second, bind_data->genotype_filter, "read_pfile");
		} else if (has_gr) {
			bool inc_missing = false;
			RangeFilter range = ParseRangeFilter(gr_it->second, "genotype_range", 0.0, 2.0, "read_pfile", &inc_missing);
			bind_data->genotype_filter.SetFromRange(range, inc_missing);
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
		// Sample orient across one OR MORE sources: pre-read every source's
		// genotypes into a single matrix concatenated along the variant axis (in
		// file/shard order). Samples are shared across sources by contract (see
		// combine_samples). All three phases below (count filter, schema, pre-read)
		// iterate sources; the matrix is addressed by GLOBAL effective-variant
		// position = variant_offsets[src] + local_ev.
		plink2::PglErr err = plink2::kPglRetSuccess;
		// Apply the count filter PER SOURCE before building schema, so the ARRAY
		// dimension reflects the filtered total. genotype_range_all_pass accumulates
		// globally across sources in file order, staying aligned with the matrix.
		vector<bool> genotype_range_all_pass;
		if (bind_data->count_filter.HasFilter() || bind_data->genotype_filter.active) {
			for (auto &src : bind_data->sources) {
				// Init temporary PgenReader for count filtering
				plink2::PgenFileInfo cf_pgfi;
				plink2::PreinitPgfi(&cf_pgfi);

				char cf_errstr[plink2::kPglErrstrBufBlen];
				plink2::PgenHeaderCtrl cf_hdr;
				uintptr_t cf_pgfi_alloc_ct = 0;

				err = plink2::PgfiInitPhase1(src.pgen_path.c_str(), nullptr, src.raw_variant_ct,
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

				err = plink2::PgrInit(src.pgen_path.c_str(), cf_max_vrec, &cf_pgfi, &cf_pgr,
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

				// Iterate this source's effective variants and filter by count
				uint32_t pre_filter_ct = src.EffectiveVariantCt();
				vector<uint32_t> filtered_indices;
				filtered_indices.reserve(pre_filter_ct);

				for (uint32_t ev = 0; ev < pre_filter_ct; ev++) {
					uint32_t vidx = src.has_effective_variant_list ? src.effective_variant_indices[ev] : ev;

					STD_ARRAY_DECL(uint32_t, 4, genocounts);
					plink2::PglErr cf_err =
					    plink2::PgrGetCounts(cf_si, cf_iv, cf_pssi, cf_sc, vidx, &cf_pgr, genocounts);
					if (cf_err != plink2::kPglRetSuccess) {
						plink2::PglErr ce = plink2::kPglRetSuccess;
						plink2::CleanupPgr(&cf_pgr, &ce);
						ce = plink2::kPglRetSuccess;
						plink2::CleanupPgfi(&cf_pgfi, &ce);
						throw IOException("read_pfile: PgrGetCounts failed for variant %u during count filter", vidx);
					}

					auto pf =
					    CheckPreDecompFilters(bind_data->count_filter, bind_data->genotype_filter, genocounts, cf_sc);
					if (pf.skip) {
						continue;
					}
					filtered_indices.push_back(vidx);
					genotype_range_all_pass.push_back(pf.all_pass);
				}

				src.effective_variant_indices = std::move(filtered_indices);
				src.has_effective_variant_list = true;

				// Cleanup
				{
					plink2::PglErr ce = plink2::kPglRetSuccess;
					plink2::CleanupPgr(&cf_pgr, &ce);
					ce = plink2::kPglRetSuccess;
					plink2::CleanupPgfi(&cf_pgfi, &ce);
				}
			}
		}

		// Rebuild variant_offsets from the POST-count-filter per-source counts, so
		// bind_data->EffectiveVariantCt() (= variant_offsets.back()) is the correct
		// concatenated total and per-source global offsets are right. Unconditional
		// (single-file included): with no count filter this reproduces the offsets
		// built at source-load time.
		{
			bind_data->variant_offsets.assign(1, 0);
			uint32_t cum = 0;
			for (auto &s : bind_data->sources) {
				cum += s.EffectiveVariantCt();
				bind_data->variant_offsets.push_back(cum);
			}
		}

		// Sample-orient mode: sample metadata columns + genotypes array/list
		auto &psam_header = bind_data->sample_metadata.header;
		for (idx_t i = 0; i < psam_header.column_names.size(); i++) {
			names.push_back(psam_header.column_names[i]);
			return_types.push_back(psam_header.column_types[i]);
			bind_data->sample_orient_col_to_psam_col.push_back(i);
		}

		// Resolve genotype output mode — dimension is the concatenated total across
		// all sources (post count filter; variant_offsets rebuilt above).
		uint32_t effective_variant_ct = bind_data->EffectiveVariantCt();
		// One variant's column/field name (ID, else CHROM:POS), for COLUMNS/STRUCT.
		auto variant_col_name = [](const PfileSource &s, uint32_t vidx) -> string {
			auto id = s.variants.GetId(vidx);
			return id.empty() ? (s.variants.GetChrom(vidx) + ":" + std::to_string(s.variants.GetPos(vidx))) : id;
		};
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

			// Check for duplicate variant IDs (variant IDs can duplicate unlike sample
			// IIDs) — across ALL sources, in global (shard) order.
			std::unordered_set<string> seen_names;
			for (auto &src : bind_data->sources) {
				uint32_t src_ct = src.EffectiveVariantCt();
				for (uint32_t ev = 0; ev < src_ct; ev++) {
					uint32_t vidx = src.has_effective_variant_list ? src.effective_variant_indices[ev] : ev;
					string col_name = variant_col_name(src, vidx);
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
			}

			bind_data->sample_orient_total_cols = names.size();
		} else if (bind_data->genotype_mode == GenotypeMode::STRUCT) {
			// STRUCT mode: single genotypes column with one field per effective variant
			child_list_t<LogicalType> struct_children;
			LogicalType field_type = bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
			                         : bind_data->include_dosages ? LogicalType::DOUBLE
			                                                      : LogicalType::TINYINT;

			std::unordered_set<string> seen_names;
			for (auto &src : bind_data->sources) {
				uint32_t src_ct = src.EffectiveVariantCt();
				for (uint32_t ev = 0; ev < src_ct; ev++) {
					uint32_t vidx = src.has_effective_variant_list ? src.effective_variant_indices[ev] : ev;
					string col_name = variant_col_name(src, vidx);
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

		// Shared decode buffers below depend only on the (shared) sample counts;
		// each source opens its own PgenReader inside the pre-read loop further down.

		// Allocate genovec buffer
		uintptr_t genovec_wc = plink2::NypCtToAlignedWordCt(bind_data->raw_sample_ct);
		AlignedBuffer genovec_buf2;
		genovec_buf2.Allocate(genovec_wc * sizeof(uintptr_t));
		std::memset(genovec_buf2.ptr, 0, genovec_wc * sizeof(uintptr_t));

		vector<int8_t> tmp_bytes(bind_data->raw_sample_ct);

		// Sample subsetting: a SINGLE subset (from the shared samples) applied to
		// every source. pssi is per-reader, so it is (re)set inside the source loop.
		SampleSubset preread_subset;
		if (bind_data->has_sample_subset) {
			preread_subset = BuildSampleSubset(bind_data->raw_sample_ct, bind_data->sample_indices);
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

		// Row-level genotype_range accumulators (sample orient): track, per output
		// sample, whether it has any in-range genotype and any missing genotype
		// across the effective variants. Computed from TRUE genotype values, before
		// any per-element null-out, so aggregate counts and the keep decision are
		// not confused by the -9 sentinel. dosages are incompatible with
		// genotype_range (rejected at parse time), so only the non-dosage paths
		// populate these.
		const bool row_filter_active = bind_data->genotype_filter.active && !bind_data->include_dosages;
		vector<uint8_t> sample_in_range;
		vector<uint8_t> sample_has_missing;
		if (row_filter_active) {
			sample_in_range.assign(output_sample_ct, 0);
			sample_has_missing.assign(output_sample_ct, 0);
		}

		// Pre-read every source into its slice of the global matrix. Source src_idx
		// occupies matrix rows [variant_offsets[src_idx], variant_offsets[src_idx+1]).
		for (idx_t src_idx = 0; src_idx < bind_data->sources.size(); src_idx++) {
			PfileSource &src = bind_data->sources[src_idx];
			uint32_t global_offset = bind_data->variant_offsets[src_idx];

			// Temporary PgenReader for THIS source.
			plink2::PgenFileInfo tmp_pgfi;
			plink2::PreinitPgfi(&tmp_pgfi);
			char errstr_buf2[plink2::kPglErrstrBufBlen];
			plink2::PgenHeaderCtrl header_ctrl2;
			uintptr_t pgfi_alloc_ct2 = 0;
			err = plink2::PgfiInitPhase1(src.pgen_path.c_str(), nullptr, src.raw_variant_ct, bind_data->raw_sample_ct,
			                             &header_ctrl2, &tmp_pgfi, &pgfi_alloc_ct2, errstr_buf2);
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
			err =
			    plink2::PgrInit(src.pgen_path.c_str(), max_vrec2, &tmp_pgfi, &tmp_pgr, pgr_alloc2.As<unsigned char>());
			if (err != plink2::kPglRetSuccess) {
				plink2::PglErr ce = plink2::kPglRetSuccess;
				plink2::CleanupPgr(&tmp_pgr, &ce);
				ce = plink2::kPglRetSuccess;
				plink2::CleanupPgfi(&tmp_pgfi, &ce);
				throw IOException("read_pfile: PgrInit failed for sample-orient pre-read");
			}
			plink2::PgrSampleSubsetIndex pssi2;
			if (bind_data->has_sample_subset) {
				plink2::PgrSetSampleSubsetIndex(preread_subset.CumulativePopcounts(), &tmp_pgr, &pssi2);
			} else {
				plink2::PgrClearSampleSubsetIndex(&tmp_pgr, &pssi2);
			}

			uint32_t src_effective_ct = src.EffectiveVariantCt();
			for (uint32_t local_ev = 0; local_ev < src_effective_ct; local_ev++) {
				uint32_t ev = global_offset + local_ev; // global matrix row
				uint32_t vidx = src.has_effective_variant_list ? src.effective_variant_indices[local_ev] : local_ev;

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
						throw IOException("read_pfile: PgrGetD failed for variant %u during sample-orient pre-read",
						                  vidx);
					}
					plink2::Dosage16ToDoublesMinus9(
					    genovec_buf2.As<uintptr_t>(), preread_dosage_present.As<uintptr_t>(),
					    preread_dosage_main.As<uint16_t>(), output_sample_ct, dosage_ct, preread_dosage_doubles.data());
					bind_data->dosage_matrix[ev].assign(preread_dosage_doubles.begin(),
					                                    preread_dosage_doubles.begin() + output_sample_ct);
				} else if (bind_data->include_phased) {
					uint32_t phasepresent_ct = 0;
					// PgrGetP leaves the phase buffers untouched for variants without an
					// hphase track; zero per-call so unphased hets decode deterministically
					// as REF|ALT [0,1] and never inherit stale phase bits. cachealigned_malloc
					// does not zero, so a one-time memset at allocation is insufficient.
					uintptr_t phase_byte_ct = plink2::BitCtToAlignedWordCt(output_sample_ct) * sizeof(uintptr_t);
					std::memset(preread_phasepresent.ptr, 0, phase_byte_ct);
					std::memset(preread_phaseinfo.ptr, 0, phase_byte_ct);
					err = plink2::PgrGetP(si_ptr, pssi2, output_sample_ct, vidx, &tmp_pgr, genovec_buf2.As<uintptr_t>(),
					                      preread_phasepresent.As<uintptr_t>(), preread_phaseinfo.As<uintptr_t>(),
					                      &phasepresent_ct);
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
					// Row-level keep accumulation (phased: diploid sum) from true values.
					if (row_filter_active) {
						for (uint32_t s = 0; s < output_sample_ct; s++) {
							int8_t a1 = preread_phased_pairs[s * 2];
							if (a1 == -9) {
								sample_has_missing[s] = 1;
							} else if (bind_data->genotype_filter.AllowsCall(
							               static_cast<double>(a1 + preread_phased_pairs[s * 2 + 1]))) {
								sample_in_range[s] = 1;
							}
						}
					}
					// Apply genotype_range per-element filtering (phased: check diploid sum)
					if (bind_data->genotype_filter.active && !genotype_range_all_pass.empty() &&
					    !genotype_range_all_pass[ev]) {
						for (uint32_t s = 0; s < output_sample_ct; s++) {
							int8_t a1 = preread_phased_pairs[s * 2];
							int8_t a2 = preread_phased_pairs[s * 2 + 1];
							if (a1 != -9 && !bind_data->genotype_filter.AllowsCall(static_cast<double>(a1 + a2))) {
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
						throw IOException("read_pfile: PgrGet failed for variant %u during sample-orient pre-read",
						                  vidx);
					}
					plink2::GenoarrToBytesMinus9(genovec_buf2.As<uintptr_t>(), output_sample_ct, tmp_bytes.data());
					// Row-level keep accumulation from true values, before any null-out.
					if (row_filter_active) {
						for (uint32_t s = 0; s < output_sample_ct; s++) {
							int8_t geno = tmp_bytes[s];
							if (geno == -9) {
								sample_has_missing[s] = 1;
							} else if (bind_data->genotype_filter.AllowsCall(static_cast<double>(geno))) {
								sample_in_range[s] = 1;
							}
						}
					}
					// Apply genotype_range per-element null-out for per-element output modes only.
					// Aggregate modes (counts/stats) must keep true genotype values so that
					// out-of-range calls are not miscounted as missing (-9) in the counts struct;
					// genotype_range acts purely as a row filter there.
					if (bind_data->genotype_filter.active && !genotype_range_all_pass.empty() &&
					    !genotype_range_all_pass[ev] && !IsAggregateGenotypeMode(bind_data->genotype_mode)) {
						for (uint32_t s = 0; s < output_sample_ct; s++) {
							int8_t geno = tmp_bytes[s];
							if (geno != -9 && !bind_data->genotype_filter.AllowsCall(static_cast<double>(geno))) {
								tmp_bytes[s] = -9;
							}
						}
					}
					bind_data->genotype_matrix[ev].assign(tmp_bytes.begin(), tmp_bytes.begin() + output_sample_ct);
				}
			}

			// Cleanup this source's temporary pgenlib state
			{
				plink2::PglErr ce = plink2::kPglRetSuccess;
				plink2::CleanupPgr(&tmp_pgr, &ce);
				ce = plink2::kPglRetSuccess;
				plink2::CleanupPgfi(&tmp_pgfi, &ce);
			}
		}

		// Build the row-level keep list: a sample is emitted iff it has at least
		// one genotype in range (or a missing genotype when include_missing is set).
		// The scan iterates this list instead of all samples, so only surviving
		// rows materialize.
		if (row_filter_active) {
			bind_data->sample_keep.reserve(output_sample_ct);
			for (uint32_t s = 0; s < output_sample_ct; s++) {
				if (sample_in_range[s] || (bind_data->genotype_filter.include_missing && sample_has_missing[s])) {
					bind_data->sample_keep.push_back(s);
				}
			}
			bind_data->has_sample_keep = true;
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

//! Projection-aware load of the parquet-source psam columns into the global state.
//! Determines which psam columns the query projects (from column_ids via the orient
//! column maps), queries ONLY those from the source, and records chunk offsets. A
//! query that projects no psam column (e.g. counts-only) builds nothing.
static void BuildProjectedPsamCdc(ClientContext &context, const PfileBindData &bind_data, PfileGlobalState &gstate,
                                  const vector<column_t> &column_ids) {
	auto &meta = bind_data.sample_metadata;
	if (!meta.from_parquet) {
		return; // text path serves values from meta.rows
	}
	if (bind_data.orient_mode != OrientMode::SAMPLE && bind_data.orient_mode != OrientMode::GENOTYPE) {
		return;
	}

	const idx_t ncol = meta.header.column_names.size();
	gstate.psam_col_to_cdc.assign(ncol, DConstants::INVALID_INDEX);

	// Collect the projected psam header columns (dedup, then sort ascending).
	vector<idx_t> projected;
	for (auto col_id : column_ids) {
		if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}
		idx_t psam_col = DConstants::INVALID_INDEX;
		if (bind_data.orient_mode == OrientMode::SAMPLE) {
			// psam metadata columns are always the leading output columns (indices
			// 0..N-1); the genotypes/per-variant columns are beyond that, in every
			// genotypes mode (array/list/counts/stats/struct/columns).
			if (col_id < bind_data.sample_orient_col_to_psam_col.size()) {
				psam_col = bind_data.sample_orient_col_to_psam_col[col_id];
			}
		} else { // GENOTYPE
			if (col_id >= bind_data.geno_sample_col_start && col_id < bind_data.geno_genotype_col) {
				idx_t rel = col_id - bind_data.geno_sample_col_start;
				if (rel < bind_data.geno_sample_col_to_psam_col.size()) {
					psam_col = bind_data.geno_sample_col_to_psam_col[rel];
				}
			}
		}
		if (psam_col != DConstants::INVALID_INDEX && psam_col < ncol &&
		    gstate.psam_col_to_cdc[psam_col] == DConstants::INVALID_INDEX) {
			gstate.psam_col_to_cdc[psam_col] = 0; // mark; real CDC index assigned below
			projected.push_back(psam_col);
		}
	}
	if (projected.empty()) {
		return; // no psam column projected — nothing to load
	}
	std::sort(projected.begin(), projected.end());
	for (idx_t i = 0; i < projected.size(); i++) {
		gstate.psam_col_to_cdc[projected[i]] = i;
	}

	auto quote = [](const string &n) {
		return "\"" + StringUtil::Replace(n, "\"", "\"\"") + "\"";
	};
	string cols;
	for (idx_t i = 0; i < projected.size(); i++) {
		cols += (i ? ", " : "") + quote(meta.header.column_names[projected[i]]);
	}
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto escaped = StringUtil::Replace(meta.source_path, "'", "''");
	auto result = conn.Query("SELECT " + cols + " FROM '" + escaped + "'");
	if (result->HasError()) {
		throw IOException("read_pfile: failed to load psam columns from '%s': %s", meta.source_path,
		                  result->GetError());
	}

	auto &collection = result->Collection();
	gstate.psam_chunk_start_row.push_back(0);
	idx_t running = 0;
	{
		ColumnDataScanState scan_state;
		collection.InitializeScan(scan_state);
		DataChunk size_chunk;
		collection.InitializeScanChunk(scan_state, size_chunk);
		while (collection.Scan(scan_state, size_chunk)) {
			running += size_chunk.size();
			gstate.psam_chunk_start_row.push_back(running);
		}
	}
	gstate.psam_cdc = result->TakeCollection();
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

// Multi-file sub-batch sizes (each batch bounded to one source). Genotype orient
// fans each variant out to N sample rows, so it uses a smaller sub-batch.
static constexpr uint32_t MULTI_VARIANT_SUBBATCH = 1000;
static constexpr uint32_t MULTI_GENOTYPE_SUBBATCH = 64;

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
		// A genotype filter must apply regardless of projection: the per-sample
		// row-skip and its genotype decode are gated on need_genotypes, so force it
		// on when the filter is active even if the genotype column is not selected.
		// Otherwise `SELECT IID FROM ... orient:='genotype', include_genotypes:=[...]`
		// (or any COUNT(*)) would emit every sample of each surviving variant unfiltered.
		if (bind_data.genotype_filter.active) {
			state->need_genotypes = true;
		}
	} else if (bind_data.orient_mode == OrientMode::SAMPLE) {
		state->total_count = bind_data.has_sample_keep ? static_cast<uint32_t>(bind_data.sample_keep.size())
		                                               : bind_data.OutputSampleCt();
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

	// Multi-file variant/genotype orient: precompute scan batches, each bounded to a
	// single source and sub-batched so no batch spans a file boundary. Threads claim
	// batches atomically (next_idx as a batch index); a thread reopens its reader at
	// most once per claimed batch. Sample orient is single-file only (guarded in bind).
	if (bind_data.sources.size() > 1 && bind_data.orient_mode != OrientMode::SAMPLE) {
		// Genotype orient fans each variant out to N sample rows → smaller sub-batches.
		const uint32_t sub =
		    (bind_data.orient_mode == OrientMode::GENOTYPE) ? MULTI_GENOTYPE_SUBBATCH : MULTI_VARIANT_SUBBATCH;
		for (uint32_t s = 0; s < bind_data.sources.size(); s++) {
			uint32_t n = bind_data.sources[s].EffectiveVariantCt();
			for (uint32_t start = 0; start < n; start += sub) {
				state->batches.push_back(ScanBatch {s, start, std::min(start + sub, n)});
			}
		}
	}

	// Projection-aware psam column load (parquet source): materialize only the psam
	// columns this query actually selects.
	BuildProjectedPsamCdc(context, bind_data, *state, input.column_ids);

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

//! Open (or reopen) this thread's PgenFileInfo/PgenReader on sources[source_idx].
//! The one-time per-thread buffers (genovec/phase/dosage/sample_include/
//! cumulative_popcounts) are sized off the shared raw_sample_ct and are NOT touched
//! here — they are built once in PfileInitLocal and reused across source swaps
//! (identical samples by contract). Only pgfi/pgr are torn down and rebuilt, and
//! pssi is re-bound to the new pgr. No-op if already open on source_idx.
static void OpenSourceReader(PfileLocalState &state, const PfileBindData &bind_data, idx_t source_idx) {
	if (state.initialized && state.current_source_idx == source_idx) {
		return;
	}
	if (state.initialized) {
		// Tear down the reader currently open on another source (pgr before pgfi).
		plink2::PglErr ce = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&state.pgr, &ce);
		plink2::CleanupPgfi(&state.pgfi, &ce);
		state.initialized = false;
	}

	const PfileSource &src = bind_data.sources[source_idx];

	plink2::PreinitPgfi(&state.pgfi);
	plink2::PreinitPgr(&state.pgr);

	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err =
	    plink2::PgfiInitPhase1(src.pgen_path.c_str(), nullptr, src.raw_variant_ct, bind_data.raw_sample_ct,
	                           &header_ctrl, &state.pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&state.pgfi, &cleanup_err);
		throw IOException("read_pfile: thread init failed (phase 1) for '%s': %s", src.pgen_path, errstr_buf);
	}

	if (pgfi_alloc_cacheline_ct > 0) {
		state.pgfi_alloc_buf.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	} else {
		state.pgfi_alloc_buf.Reset();
	}

	uint32_t max_vrec_width = 0;
	uintptr_t pgr_alloc_cacheline_ct = 0;
	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, state.pgfi.raw_variant_ct, &max_vrec_width, &state.pgfi,
	                             state.pgfi_alloc_buf.As<unsigned char>(), &pgr_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&state.pgfi, &cleanup_err);
		throw IOException("read_pfile: thread init failed (phase 2) for '%s': %s", src.pgen_path, errstr_buf);
	}

	if (pgr_alloc_cacheline_ct > 0) {
		state.pgr_alloc_buf.Allocate(pgr_alloc_cacheline_ct * plink2::kCacheline);
	} else {
		state.pgr_alloc_buf.Reset();
	}

	err = plink2::PgrInit(src.pgen_path.c_str(), max_vrec_width, &state.pgfi, &state.pgr,
	                      state.pgr_alloc_buf.As<unsigned char>());
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&state.pgr, &cleanup_err);
		cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&state.pgfi, &cleanup_err);
		throw IOException("read_pfile: PgrInit failed for '%s'", src.pgen_path);
	}

	// Re-bind the sample-subset index to the (new) pgr. cumulative_popcounts_buf was
	// filled once in PfileInitLocal and is identical across sources.
	if (bind_data.has_sample_subset) {
		plink2::PgrSetSampleSubsetIndex(state.cumulative_popcounts_buf.As<uint32_t>(), &state.pgr, &state.pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state.pgr, &state.pssi);
	}

	state.current_source_idx = source_idx;
	state.initialized = true;
}

static unique_ptr<LocalTableFunctionState> PfileInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                          GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PfileBindData>();
	auto &gstate = global_state->Cast<PfileGlobalState>();
	auto state = make_uniq<PfileLocalState>();

	if (!gstate.need_pgen_reader || bind_data.orient_mode == OrientMode::SAMPLE) {
		// Sample-orient mode uses pre-read genotype matrix — no per-thread PgenReader needed
		return std::move(state);
	}

	// --- Allocate one-time per-thread buffers (sized off shared raw_sample_ct) ---
	// These are reused across all sources (identical samples by contract); only the
	// pgfi/pgr reader swaps at a file boundary.
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

	// Build the sample-subset bit vector + cumulative popcounts once (shared samples).
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
	}

	// Open the first source (single-file: the only one). Multi-file scans reopen on a
	// source boundary via OpenSourceReader.
	OpenSourceReader(*state, bind_data, 0);

	return std::move(state);
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

	// Emit one variant's row at output position `rows_emitted` (captured), reading
	// genotypes from `source` (the caller must have this source's reader open when
	// gstate.need_pgen_reader). Returns false if the variant was skipped by a count/
	// genotype pre-decompression filter (no row emitted). Shared by the single-file
	// and multi-file claim loops below.
	auto emit_variant_row = [&](const PfileSource &source, uint32_t local_pos) -> bool {
		uint32_t vidx = source.ResolveVariantIdx(local_pos);

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
				return false;
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
				plink2::Dosage16ToDoublesMinus9(
				    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
				    lstate.dosage_main_buf.As<uint16_t>(), output_sample_ct, dosage_ct, lstate.dosage_doubles.data());
			} else if (bind_data.include_phased) {
				// PgrGetP does NOT touch the phase buffers for variants without an
				// hphase track (it early-returns after setting phasepresent_ct=0).
				// Zero them per-call so unphased hets decode deterministically as
				// REF|ALT [0,1] and never inherit stale phase bits from a prior
				// (phased) variant. cachealigned_malloc does not zero, so a
				// one-time memset at allocation is insufficient.
				uintptr_t phase_byte_ct = plink2::BitCtToAlignedWordCt(output_sample_ct) * sizeof(uintptr_t);
				std::memset(lstate.phasepresent_buf.ptr, 0, phase_byte_ct);
				std::memset(lstate.phaseinfo_buf.ptr, 0, phase_byte_ct);
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
				auto val = source.variants.GetChrom(vidx);
				FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				break;
			}
			case PfileBindData::POS_COL: {
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = source.variants.GetPos(vidx);
				break;
			}
			case PfileBindData::ID_COL: {
				auto val = source.variants.GetId(vidx);
				if (val.empty()) {
					FlatVector::SetNull(vec, rows_emitted, true);
				} else {
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				}
				break;
			}
			case PfileBindData::REF_COL: {
				auto val = source.variants.GetRef(vidx);
				FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				break;
			}
			case PfileBindData::ALT_COL: {
				auto val = source.variants.GetAlt(vidx);
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
							if (geno == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
							                   !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno)))) {
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
							if (geno == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
							                   !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno)))) {
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
					if (!bind_data.genotype_filter.active) {
						// No genotype range filter: use shared output utility
						FillGenotypeVector(vec, rows_emitted, bind_data.genotype_mode, output_sample_ct, nullptr,
						                   lstate.phased_pairs.data(), true);
					} else if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
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
							if (a1 == -9 || (!geno_range_all_pass &&
							                 !bind_data.genotype_filter.AllowsCall(static_cast<double>(a1 + a2)))) {
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
							if (a1 == -9 || (!geno_range_all_pass &&
							                 !bind_data.genotype_filter.AllowsCall(static_cast<double>(a1 + a2)))) {
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
				} else if (!bind_data.genotype_filter.active) {
					// No genotype range filter: use shared output utility
					FillGenotypeVector(vec, rows_emitted, bind_data.genotype_mode, output_sample_ct,
					                   lstate.genotype_bytes.data(), nullptr, false);
				} else {
					if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
						auto array_size = static_cast<idx_t>(output_sample_ct);
						auto &child = ArrayVector::GetEntry(vec);
						auto *child_data = FlatVector::GetData<int8_t>(child);
						auto &child_validity = FlatVector::Validity(child);

						idx_t base = rows_emitted * array_size;
						for (idx_t s = 0; s < array_size; s++) {
							int8_t geno = lstate.genotype_bytes[s];
							if (geno == -9 || (!geno_range_all_pass &&
							                   !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno)))) {
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
							if (geno == -9 || (!geno_range_all_pass &&
							                   !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno)))) {
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
							if (geno == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
							                   !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno)))) {
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

		return true;
	};

	if (bind_data.sources.size() == 1) {
		// Single-file: original dynamic fetch_add claim loop over sources[0] (byte-for-byte
		// identical to the pre-multi-file behavior).
		const PfileSource &source = bind_data.Primary();
		while (rows_emitted < STANDARD_VECTOR_SIZE) {
			uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
			uint32_t claim_size = std::min(PFILE_BATCH_SIZE, remaining_capacity);
			uint32_t batch_start = gstate.next_idx.fetch_add(claim_size);
			if (batch_start >= total_variants) {
				break;
			}
			uint32_t batch_end = std::min(batch_start + claim_size, total_variants);
			for (uint32_t effective_pos = batch_start; effective_pos < batch_end; effective_pos++) {
				if (emit_variant_row(source, effective_pos)) {
					rows_emitted++;
				}
			}
		}
	} else {
		// Multi-file: claim precomputed batches (each bounded to one source), reopening the
		// per-thread reader only at a source boundary (OpenSourceReader is a no-op otherwise).
		while (rows_emitted < STANDARD_VECTOR_SIZE) {
			if (lstate.mf_need_claim) {
				lstate.mf_batch = gstate.next_idx.fetch_add(1);
				if (lstate.mf_batch >= gstate.batches.size()) {
					break;
				}
				lstate.mf_local = gstate.batches[lstate.mf_batch].local_start;
				lstate.mf_need_claim = false;
			}
			const ScanBatch &batch = gstate.batches[lstate.mf_batch];
			const PfileSource &source = bind_data.sources[batch.source_idx];
			if (gstate.need_pgen_reader) {
				OpenSourceReader(lstate, bind_data, batch.source_idx);
			}
			while (lstate.mf_local < batch.local_end && rows_emitted < STANDARD_VECTOR_SIZE) {
				if (emit_variant_row(source, lstate.mf_local)) {
					rows_emitted++;
				}
				lstate.mf_local++;
			}
			if (lstate.mf_local >= batch.local_end) {
				lstate.mf_need_claim = true;
			}
		}
	}

	CompatSetOutputCardinality(output, rows_emitted);
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

//! Emit one psam column value for a sample into an output vector:
//!  - text path: the std::string in PfileSampleMetadata::rows.
//!  - parquet path: read from the projection-aware CDC in the global state via a
//!    per-thread cached FetchChunk. The CDC holds ONLY the projected psam columns
//!    (psam_col_to_cdc maps header col -> CDC col), so a wide psam never fetches
//!    unused columns. VARCHAR columns read string_t directly (no Value boxing);
//!    rarely-projected non-VARCHAR columns box the single value. All values go
//!    through FillSampleMetadataValue for identical SEX/PAT/MAT/missing handling.
static void EmitSampleMetadataColumn(const PfileBindData &bind_data, PfileGlobalState &gstate, PfileLocalState &lstate,
                                     Vector &out_vec, idx_t out_row, idx_t sample_file_idx, idx_t psam_col_idx) {
	auto &meta = bind_data.sample_metadata;

	if (!meta.from_parquet) {
		if (sample_file_idx < meta.rows.size() && psam_col_idx < meta.rows[sample_file_idx].size()) {
			FillSampleMetadataValue(out_vec, out_row, meta.rows[sample_file_idx][psam_col_idx], psam_col_idx, meta);
		} else {
			FlatVector::SetNull(out_vec, out_row, true);
		}
		return;
	}

	// Parquet path. Map the header column to its column in the projected CDC.
	if (!gstate.psam_cdc || psam_col_idx >= gstate.psam_col_to_cdc.size() ||
	    gstate.psam_col_to_cdc[psam_col_idx] == DConstants::INVALID_INDEX) {
		FlatVector::SetNull(out_vec, out_row, true);
		return;
	}
	idx_t cdc_col = gstate.psam_col_to_cdc[psam_col_idx];
	const auto &starts = gstate.psam_chunk_start_row;
	if (starts.size() < 2 || sample_file_idx >= starts.back()) {
		FlatVector::SetNull(out_vec, out_row, true);
		return;
	}
	idx_t chunk_idx = (idx_t)(std::upper_bound(starts.begin(), starts.end(), sample_file_idx) - starts.begin()) - 1;
	if (chunk_idx != lstate.psam_chunk_idx) {
		if (lstate.psam_chunk.ColumnCount() == 0) {
			gstate.psam_cdc->InitializeScanChunk(lstate.psam_chunk);
		}
		gstate.psam_cdc->FetchChunk(chunk_idx, lstate.psam_chunk);
		lstate.psam_chunk.Flatten(); // FlatVector access below; only the projected cols
		lstate.psam_chunk_idx = chunk_idx;
	}
	idx_t local = sample_file_idx - starts[chunk_idx];
	Vector &src = lstate.psam_chunk.data[cdc_col];

	// SEX/PAT/MAT and any VARCHAR-typed output column go through the string-based
	// FillSampleMetadataValue (SEX 0->NULL, PAT/MAT "0"->NULL, missing sentinels).
	bool is_special = (psam_col_idx == meta.sex_col_idx);
	for (auto p : meta.parent_col_indices) {
		if (p == psam_col_idx) {
			is_special = true;
			break;
		}
	}
	if (is_special || out_vec.GetType().id() == LogicalTypeId::VARCHAR) {
		if (src.GetType().id() == LogicalTypeId::VARCHAR) {
			if (!FlatVector::Validity(src).RowIsValid(local)) {
				FillSampleMetadataValue(out_vec, out_row, string(), psam_col_idx, meta);
			} else {
				FillSampleMetadataValue(out_vec, out_row, FlatVector::GetData<string_t>(src)[local].GetString(),
				                        psam_col_idx, meta);
			}
		} else {
			Value v = src.GetValue(local);
			FillSampleMetadataValue(out_vec, out_row, v.IsNull() ? string() : v.ToString(), psam_col_idx, meta);
		}
	} else {
		// Typed covariate column (INTEGER/DOUBLE/...): output type matches the source
		// parquet type, so copy the value natively rather than through a string.
		out_vec.SetValue(out_row, src.GetValue(local));
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

	// Decode one variant's genotypes into the per-thread buffers (reader must be open
	// on `source`). Returns false if the variant is skipped by a count/genotype
	// pre-decompression filter (no rows emitted for it). Sets lstate.geno_range_all_pass.
	auto load_variant = [&](const PfileSource &source, uint32_t vidx) -> bool {
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
				return false;
			}
			lstate.geno_range_all_pass = pf.all_pass;
		}

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
				plink2::Dosage16ToDoublesMinus9(
				    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
				    lstate.dosage_main_buf.As<uint16_t>(), output_sample_ct, dosage_ct, lstate.dosage_doubles.data());
			} else if (bind_data.include_phased) {
				// PgrGetP does NOT touch the phase buffers for variants without an
				// hphase track (it early-returns after setting phasepresent_ct=0).
				// Zero them per-call so unphased hets decode deterministically as
				// REF|ALT [0,1] and never inherit stale phase bits from a prior
				// (phased) variant. cachealigned_malloc does not zero, so a
				// one-time memset at allocation is insufficient.
				uintptr_t phase_byte_ct = plink2::BitCtToAlignedWordCt(output_sample_ct) * sizeof(uintptr_t);
				std::memset(lstate.phasepresent_buf.ptr, 0, phase_byte_ct);
				std::memset(lstate.phaseinfo_buf.ptr, 0, phase_byte_ct);
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
		}
		return true;
	};

	// Emit one (variant, sample) row at output position `rows_emitted` (captured).
	// Returns false if a genotype_range row filter skips this sample (no row emitted).
	// Assumes lstate.batch_variant_loaded reflects whether genotypes were decoded.
	auto emit_sample_row = [&](const PfileSource &source, uint32_t vidx, uint32_t sample_in_variant) -> bool {
		// When genotype_range is active, skip rows for missing/out-of-range genotypes
		if (bind_data.genotype_filter.active && gstate.need_genotypes && lstate.batch_variant_loaded) {
			if (bind_data.include_phased) {
				int8_t a1 = lstate.phased_pairs[sample_in_variant * 2];
				if (a1 == -9) {
					if (!bind_data.genotype_filter.include_missing) {
						return false;
					}
				} else if (!lstate.geno_range_all_pass) {
					int8_t a2 = lstate.phased_pairs[sample_in_variant * 2 + 1];
					if (!bind_data.genotype_filter.AllowsCall(static_cast<double>(a1 + a2))) {
						return false;
					}
				}
			} else if (!bind_data.include_dosages) {
				int8_t geno = lstate.genotype_bytes[sample_in_variant];
				if (geno == -9) {
					if (!bind_data.genotype_filter.include_missing) {
						return false;
					}
				} else if (!lstate.geno_range_all_pass &&
				           !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno))) {
					return false;
				}
			}
		}

		// Map scan-order sample index to file-order sample index.
		// When subsetting, sample_indices is sorted to match pgenlib's
		// ascending-bit-order output, so genotype_bytes[i] corresponds
		// to sample_indices[i].
		uint32_t sample_file_idx =
		    bind_data.has_sample_subset ? bind_data.sample_indices[sample_in_variant] : sample_in_variant;

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
					auto val = source.variants.GetChrom(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case PfileBindData::POS_COL: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = source.variants.GetPos(vidx);
					break;
				}
				case PfileBindData::ID_COL: {
					auto val = source.variants.GetId(vidx);
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case PfileBindData::REF_COL: {
					auto val = source.variants.GetRef(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case PfileBindData::ALT_COL: {
					auto val = source.variants.GetAlt(vidx);
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
						int8_t a1 = lstate.phased_pairs[sample_in_variant * 2];
						if (a1 == -9) {
							FlatVector::SetNull(vec, rows_emitted, true);
							allele_data[allele_base] = 0;
							allele_data[allele_base + 1] = 0;
						} else {
							allele_data[allele_base] = a1;
							allele_data[allele_base + 1] = lstate.phased_pairs[sample_in_variant * 2 + 1];
						}
					} else if (bind_data.include_dosages) {
						// Scalar DOUBLE
						double dosage = lstate.dosage_doubles[sample_in_variant];
						if (dosage == -9.0) {
							FlatVector::SetNull(vec, rows_emitted, true);
						} else {
							FlatVector::GetData<double>(vec)[rows_emitted] = dosage;
						}
					} else {
						// Scalar TINYINT
						int8_t geno = lstate.genotype_bytes[sample_in_variant];
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

				EmitSampleMetadataColumn(bind_data, gstate, lstate, vec, rows_emitted, sample_file_idx, psam_col_idx);
			}
		}
		return true;
	};

	if (bind_data.sources.size() == 1) {
		// Single-file: original 64-variant fetch_add batch claim over sources[0].
		const PfileSource &source = bind_data.Primary();
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

			uint32_t vidx = source.ResolveVariantIdx(lstate.current_variant_in_batch);

			if (!lstate.batch_variant_loaded) {
				if (!load_variant(source, vidx)) {
					lstate.current_variant_in_batch++;
					lstate.current_sample_in_variant = 0;
					continue;
				}
				lstate.batch_variant_loaded = true;
			}

			// Emit rows for samples within current variant
			while (lstate.current_sample_in_variant < output_sample_ct && rows_emitted < STANDARD_VECTOR_SIZE) {
				if (emit_sample_row(source, vidx, lstate.current_sample_in_variant)) {
					rows_emitted++;
				}
				lstate.current_sample_in_variant++;
			}

			// Advance to next variant if all samples emitted
			if (lstate.current_sample_in_variant >= output_sample_ct) {
				lstate.current_variant_in_batch++;
				lstate.current_sample_in_variant = 0;
				lstate.batch_variant_loaded = false;
			}
		}
	} else {
		// Multi-file: claim precomputed batches (each bounded to one source), reopening the
		// per-thread reader only at a source boundary. mf_local is the current variant's
		// effective position within the batch; current_sample_in_variant persists mid-variant.
		while (rows_emitted < STANDARD_VECTOR_SIZE) {
			if (lstate.mf_need_claim) {
				lstate.mf_batch = gstate.next_idx.fetch_add(1);
				if (lstate.mf_batch >= gstate.batches.size()) {
					break;
				}
				lstate.mf_local = gstate.batches[lstate.mf_batch].local_start;
				lstate.current_sample_in_variant = 0;
				lstate.batch_variant_loaded = false;
				lstate.mf_need_claim = false;
			}
			const ScanBatch &batch = gstate.batches[lstate.mf_batch];
			const PfileSource &source = bind_data.sources[batch.source_idx];
			if (gstate.need_pgen_reader) {
				OpenSourceReader(lstate, bind_data, batch.source_idx);
			}

			while (lstate.mf_local < batch.local_end && rows_emitted < STANDARD_VECTOR_SIZE) {
				uint32_t vidx = source.ResolveVariantIdx(lstate.mf_local);

				if (!lstate.batch_variant_loaded) {
					if (!load_variant(source, vidx)) {
						lstate.mf_local++;
						lstate.current_sample_in_variant = 0;
						continue;
					}
					lstate.batch_variant_loaded = true;
				}

				while (lstate.current_sample_in_variant < output_sample_ct && rows_emitted < STANDARD_VECTOR_SIZE) {
					if (emit_sample_row(source, vidx, lstate.current_sample_in_variant)) {
						rows_emitted++;
					}
					lstate.current_sample_in_variant++;
				}

				if (lstate.current_sample_in_variant >= output_sample_ct) {
					lstate.mf_local++;
					lstate.current_sample_in_variant = 0;
					lstate.batch_variant_loaded = false;
				}
			}

			if (lstate.mf_local >= batch.local_end) {
				lstate.mf_need_claim = true;
			}
		}
	}

	CompatSetOutputCardinality(output, rows_emitted);
}

// ---------------------------------------------------------------------------
// Scan: Sample-orient mode (one row per sample)
// ---------------------------------------------------------------------------

static void PfileSampleOrientScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PfileBindData>();
	auto &gstate = data_p.global_state->Cast<PfileGlobalState>();
	auto &lstate = data_p.local_state->Cast<PfileLocalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t total_samples = gstate.total_count;
	// Sample orient concatenates all sources along the variant axis. PfileBind
	// rebuilds variant_offsets AFTER the per-source count filter, so
	// bind_data.EffectiveVariantCt() (= variant_offsets.back()) is the correct
	// post-filter total matrix width across every shard.
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

		for (uint32_t claim_idx = batch_start; claim_idx < batch_end; claim_idx++) {
			// When a genotype_range row filter is active, the claimed index addresses
			// the surviving-sample list; otherwise it is the output sample position directly.
			uint32_t sample_pos = bind_data.has_sample_keep ? bind_data.sample_keep[claim_idx] : claim_idx;
			// Map output sample position to file-order sample index
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

					EmitSampleMetadataColumn(bind_data, gstate, lstate, vec, rows_emitted, sample_file_idx,
					                         psam_col_idx);
				}
			}

			rows_emitted++;
		}
	}

	CompatSetOutputCardinality(output, rows_emitted);
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
	// Two overloads: a single prefix (VARCHAR) and a list of prefixes (LIST(VARCHAR),
	// row-concatenated across variant-sharded pfiles). Both dispatch to the same
	// callbacks; bind detects the list via ResolvePathList.
	auto add_named_params = [](TableFunction &fn) {
		fn.projection_pushdown = true;
		fn.named_parameters["pgen"] = LogicalType::VARCHAR;
		fn.named_parameters["pvar"] = LogicalType::VARCHAR;
		fn.named_parameters["psam"] = LogicalType::VARCHAR;
		fn.named_parameters["orient"] = LogicalType::VARCHAR;
		fn.named_parameters["dosages"] = LogicalType::BOOLEAN;
		fn.named_parameters["phased"] = LogicalType::BOOLEAN;
		fn.named_parameters["region"] = LogicalType::VARCHAR;
		// Accept ANY for samples and variants — type dispatch handled in bind
		fn.named_parameters["samples"] = LogicalType::ANY;
		fn.named_parameters["variants"] = LogicalType::ANY;
		fn.named_parameters["genotypes"] = LogicalType::VARCHAR;
		fn.named_parameters["af_range"] = LogicalType::ANY;
		fn.named_parameters["ac_range"] = LogicalType::ANY;
		fn.named_parameters["genotype_range"] = LogicalType::ANY;
		fn.named_parameters["include_genotypes"] = LogicalType::LIST(LogicalType::VARCHAR);
		fn.named_parameters["combine_samples"] = LogicalType::VARCHAR;
	};

	TableFunctionSet set("read_pfile");
	TableFunction one("read_pfile", {LogicalType::VARCHAR}, PfileScan, PfileBind, PfileInitGlobal, PfileInitLocal);
	add_named_params(one);
	set.AddFunction(one);
	TableFunction many("read_pfile", {LogicalType::LIST(LogicalType::VARCHAR)}, PfileScan, PfileBind, PfileInitGlobal,
	                   PfileInitLocal);
	add_named_params(many);
	set.AddFunction(many);
	loader.RegisterFunction(set);
}

} // namespace duckdb
