#include "pgen_reader.hpp"
#include "duckdb_compat.hpp"
#include "plink_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/function/function_set.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <unordered_set>

namespace duckdb {

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

//! One variant-sharded pgen source. A single-file read has exactly one; a
//! LIST(VARCHAR) read has one per pgen path, row-concatenated in list order. All
//! sources share identical samples by contract (same .psam) — only the variants
//! differ. read_pgen is variant-orient only, so this is pure row-concat. Sample
//! metadata lives on PgenBindData (from sources[0]).
struct PgenSource {
	string pgen_path;
	string pvar_path;

	// Variant metadata (offset-indexed for on-demand parsing)
	VariantMetadataIndex variants;

	// pgenlib header: this file's variant count. raw_sample_ct is shared (PgenBindData).
	uint32_t raw_variant_ct = 0;

	// Effective variant list for this source (from the `variants` filter, single-file
	// only). When has_effective_variant_list is false, scan all raw_variant_ct.
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

struct PgenBindData : public TableFunctionData {
	// One or more variant-sharded sources, row-concatenated in list order.
	// Single-file reads have exactly one; sources[0] is authoritative for the shared
	// sample metadata (identical .psam across shards by contract).
	vector<PgenSource> sources;

	// Cumulative effective-variant offsets across sources (size sources.size()+1).
	// Global effective position g in [variant_offsets[i], variant_offsets[i+1]) maps
	// to source i, local position g - variant_offsets[i]. Used by multi-file scan.
	vector<uint32_t> variant_offsets;

	// Shared sample metadata (from sources[0]; .psam is optional for read_pgen).
	string psam_path;
	SampleInfo sample_info;
	bool has_sample_info = false;
	uint32_t sample_ct = 0; // from .pgen header if no .psam

	// pgenlib header: sample count (shared — identical across all sources by contract)
	uint32_t raw_sample_ct = 0;

	// Column layout for projection pushdown
	// Fixed metadata: CHROM(0), POS(1), ID(2), REF(3), ALT(4)
	static constexpr idx_t GENOTYPES_COL_IDX = 5;

	// Options
	bool include_dosages = false;
	bool include_phased = false;
	GenotypeMode genotype_mode = GenotypeMode::ARRAY;

	// Sample subsetting
	bool has_sample_subset = false;
	vector<uint32_t> sample_indices; // 0-based indices into .pgen sample order
	uint32_t subset_sample_ct = 0;

	// Count-based filtering (af_range, ac_range)
	CountFilter count_filter;
	unique_ptr<SampleSubset> count_filter_subset;

	// Genotype range filtering (genotype_range)
	GenotypeRangeFilter genotype_filter;

	// Variant filtering (single-file only — resolved into the source's effective list)
	bool has_variant_filter = false;

	// Columns mode layout (genotypes := 'columns')
	vector<string> genotype_column_names;     // IIDs for column names
	idx_t columns_mode_first_geno_col = 0;    // first genotype column index in schema
	uint32_t columns_mode_geno_col_count = 0; // number of genotype columns

	//! Total effective variants across all sources (global scan count).
	uint32_t EffectiveVariantCt() const {
		return variant_offsets.empty() ? 0 : variant_offsets.back();
	}

	//! Convenience accessor for the primary (and, single-file, only) source.
	PgenSource &Primary() {
		return sources[0];
	}
	const PgenSource &Primary() const {
		return sources[0];
	}
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

//! A unit of scan work bounded to a single source (multi-file). local_start/local_end
//! are 0-based effective positions within that source. Threads claim batches atomically
//! (next_variant_idx as a batch index); a batch never spans a file, so a thread reopens
//! the pgen reader at most once per claimed batch.
struct ScanBatch {
	uint32_t source_idx;
	uint32_t local_start;
	uint32_t local_end;
};

struct PgenGlobalState : public GlobalTableFunctionState {
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t total_variants = 0;

	// Multi-file (sources.size() > 1): precomputed batches, each bounded to one source.
	// Claimed via next_variant_idx (as a batch index). Empty for single-file reads, which
	// keep the original dynamic fetch_add claim loop.
	vector<ScanBatch> batches;

	// Projection flags (computed once in init)
	bool need_genotypes = false;
	bool need_pgen_reader = false;
	vector<column_t> column_ids;
	uint32_t max_threads_config = 0;

	idx_t MaxThreads() const override {
		return ApplyMaxThreadsCap(total_variants / 1000 + 1, max_threads_config);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PgenLocalState : public LocalTableFunctionState {
	// Per-thread PgenFileInfo — must outlive the PgenReader since pgr
	// holds a reference to pgfi's shared state.
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;
	AlignedBuffer genovec_buf;
	AlignedBuffer sample_include_buf;
	AlignedBuffer cumulative_popcounts_buf;

	// Expanded genotype buffer (int8, one per sample, -9 = missing)
	vector<int8_t> genotype_bytes;

	// Phase buffers (for phased := true)
	AlignedBuffer phasepresent_buf;
	AlignedBuffer phaseinfo_buf;
	vector<int8_t> phased_pairs; // sample_ct * 2 entries
	uint32_t phasepresent_ct = 0;

	// Dosage buffers (for dosages := true)
	AlignedBuffer dosage_present_buf;
	AlignedBuffer dosage_main_buf;
	vector<double> dosage_doubles;

	// Sample subset index (for PgrGet)
	plink2::PgrSampleSubsetIndex pssi;

	bool initialized = false;

	// Multi-file: index of the source this thread's pgr/pgfi are currently open on
	// (DConstants::INVALID_INDEX = none open yet). Reopened on a source boundary.
	idx_t current_source_idx = DConstants::INVALID_INDEX;

	// Multi-file batch-claim position. mf_batch is the claimed index into
	// PgenGlobalState::batches; mf_local is the next effective position within that
	// batch's [local_start, local_end). mf_need_claim triggers claiming a fresh batch.
	uint32_t mf_batch = 0;
	uint32_t mf_local = 0;
	bool mf_need_claim = true;

	~PgenLocalState() {
		if (initialized) {
			// PgenReader must be cleaned up before PgenFileInfo
			plink2::PglErr reterr = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr, &reterr);
			plink2::CleanupPgfi(&pgfi, &reterr);
			// AlignedBuffer destructors handle the aligned allocs
		}
	}
};

// ---------------------------------------------------------------------------
// Per-source loading (companion discovery + pgenlib header + variant metadata)
// ---------------------------------------------------------------------------

//! Discover the .pvar companion, read the .pgen header, and load variant metadata for
//! ONE source. `override_pvar` supplies an explicit path (single-file only; empty for
//! shards in a list). Outputs this file's sample count via `raw_sample_ct_out` for the
//! caller's cross-source equality check.
static PgenSource LoadPgenSource(ClientContext &context, FileSystem &fs, const string &pgen_path,
                                 const string &override_pvar, uint32_t &raw_sample_ct_out) {
	PgenSource src;
	src.pgen_path = pgen_path;
	src.pvar_path = override_pvar;

	// --- Auto-discover .pvar/.bim companion ---
	if (src.pvar_path.empty()) {
		src.pvar_path = FindCompanionFileWithParquet(context, fs, src.pgen_path, {".pvar", ".bim"});
		if (src.pvar_path.empty()) {
			throw InvalidInputException("read_pgen: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            src.pgen_path);
		}
	}

	// --- Read the .pgen header (counts) ---
	plink2::PgenFileInfo pgfi;
	plink2::PreinitPgfi(&pgfi);

	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err = plink2::PgfiInitPhase1(src.pgen_path.c_str(), nullptr, // no .pgi file
	                                            UINT32_MAX,                     // infer raw_variant_ct
	                                            UINT32_MAX,                     // infer raw_sample_ct
	                                            &header_ctrl, &pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("read_pgen: failed to open '%s': %s", src.pgen_path, errstr_buf);
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

	// Clean up pgfi — we only needed it to get counts and alloc sizes.
	plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
	plink2::CleanupPgfi(&pgfi, &cleanup_err);

	if (err != plink2::kPglRetSuccess) {
		throw IOException("read_pgen: failed to initialize '%s' (phase 2): %s", src.pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	src.variants = LoadVariantMetadata(context, src.pvar_path, "read_pgen");

	if (src.variants.variant_ct != src.raw_variant_ct) {
		throw InvalidInputException("read_pgen: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            src.raw_variant_ct, src.pvar_path,
		                            static_cast<unsigned long long>(src.variants.variant_ct));
	}

	return src;
}

// ---------------------------------------------------------------------------
// Init local helper: open (or reopen) this thread's reader on a given source
// ---------------------------------------------------------------------------

//! Open (or reopen) this thread's PgenFileInfo/PgenReader on sources[source_idx].
//! The one-time per-thread buffers (genovec/phase/dosage/sample_include/
//! cumulative_popcounts) are sized off the shared raw_sample_ct and are NOT touched
//! here — they are built once in PgenInitLocal and reused across source swaps
//! (identical samples by contract). Only pgfi/pgr are torn down and rebuilt, and pssi
//! is re-bound to the new pgr. No-op if already open on source_idx.
static void OpenSourceReader(PgenLocalState &state, const PgenBindData &bind_data, idx_t source_idx) {
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

	const PgenSource &src = bind_data.sources[source_idx];

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
		throw IOException("read_pgen: thread init failed (phase 1) for '%s': %s", src.pgen_path, errstr_buf);
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
		throw IOException("read_pgen: thread init failed (phase 2) for '%s': %s", src.pgen_path, errstr_buf);
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
		throw IOException("read_pgen: PgrInit failed for '%s'", src.pgen_path);
	}

	// Re-bind the sample-subset index to the (new) pgr. cumulative_popcounts_buf was
	// filled once in PgenInitLocal and is identical across sources.
	if (bind_data.has_sample_subset) {
		plink2::PgrSetSampleSubsetIndex(state.cumulative_popcounts_buf.As<uint32_t>(), &state.pgr, &state.pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state.pgr, &state.pssi);
	}

	state.current_source_idx = source_idx;
	state.initialized = true;
}

// ---------------------------------------------------------------------------
// Bind function
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PgenBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PgenBindData>();

	// First positional argument is a .pgen path (VARCHAR) or a list of .pgen paths
	// (LIST(VARCHAR), row-concatenated in list order).
	vector<string> pgen_paths = ResolvePathList(input.inputs[0], "read_pgen");
	const bool multi_file = pgen_paths.size() > 1;

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters ---
	string override_pvar, override_psam;
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			override_pvar = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			override_psam = kv.second.GetValue<string>();
		} else if (kv.first == "dosages") {
			bind_data->include_dosages = kv.second.GetValue<bool>();
		} else if (kv.first == "phased") {
			bind_data->include_phased = kv.second.GetValue<bool>();
		} else if (kv.first == "orient") {
			auto orient_val = StringUtil::Lower(kv.second.GetValue<string>());
			if (orient_val != "variant") {
				throw InvalidInputException("read_pgen: orient := '%s' is not supported "
				                            "(read_pgen only supports orient := 'variant'; "
				                            "use read_pfile for orient := 'genotype' or 'sample')",
				                            kv.second.GetValue<string>());
			}
		} else if (kv.first == "samples") {
			// Handled after pgenlib init (need raw_sample_ct)
		} else if (kv.first == "genotypes") {
			// Handled after sample count is known
		} else if (kv.first == "variants") {
			// Handled after variant metadata is loaded
		}
		// af_range, ac_range handled after pgenlib init
	}

	if (bind_data->include_dosages && bind_data->include_phased) {
		throw InvalidInputException("read_pgen: dosages and phased cannot both be true");
	}

	// --- Multi-file guards (throw before loading anything) ---
	if (multi_file) {
		if (!override_pvar.empty() || !override_psam.empty()) {
			throw InvalidInputException("read_pgen: pvar/psam overrides cannot be combined with a multi-file "
			                            "list (they can't disambiguate per shard)");
		}
		// `variants` resolves against a single file's variant index; across a multi-file
		// list neither by-ID nor by-index resolves correctly yet (see read_pfile for the
		// global-resolution follow-up). Guard it — use per-file reads for now.
		if (input.named_parameters.find("variants") != input.named_parameters.end()) {
			throw InvalidInputException("read_pgen: variants := [...] with a multi-file list is not yet supported");
		}
	}

	// --- Build sources (one per pgen path) ---
	for (idx_t si = 0; si < pgen_paths.size(); si++) {
		uint32_t src_sample_ct = 0;
		PgenSource src = LoadPgenSource(context, fs, pgen_paths[si], si == 0 ? override_pvar : string(), src_sample_ct);

		// Free safety check: every shard's .pgen header reports its sample_ct. All shards
		// share identical samples by contract; a genuine mismatch would silently overflow
		// the per-thread sample buffers — reject it clearly.
		if (si == 0) {
			bind_data->raw_sample_ct = src_sample_ct;
		} else if (src_sample_ct != bind_data->raw_sample_ct) {
			throw InvalidInputException(
			    "read_pgen: sample count mismatch across files: '%s' has %u samples, but expected %u "
			    "(from '%s'). All listed pgen files must share the same samples.",
			    src.pgen_path, src_sample_ct, bind_data->raw_sample_ct, bind_data->sources[0].pgen_path);
		}

		bind_data->sources.push_back(std::move(src));
	}

	// --- Discover + load sample info (from sources[0] only — identical by contract) ---
	bind_data->psam_path = override_psam;
	if (bind_data->psam_path.empty()) {
		bind_data->psam_path =
		    FindCompanionFileWithParquet(context, fs, bind_data->Primary().pgen_path, {".psam", ".fam"});
		// .psam is optional for read_pgen — if not found, we operate in index-only mode
	}

	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleMetadata(context, bind_data->psam_path);
		bind_data->has_sample_info = true;
		bind_data->sample_ct = static_cast<uint32_t>(bind_data->sample_info.sample_ct);

		if (bind_data->sample_ct != bind_data->raw_sample_ct) {
			throw InvalidInputException("read_pgen: sample count mismatch: .pgen has %u samples, "
			                            ".psam/.fam '%s' has %u samples",
			                            bind_data->raw_sample_ct, bind_data->psam_path, bind_data->sample_ct);
		}
	} else {
		bind_data->sample_ct = bind_data->raw_sample_ct;
	}

	// --- Process samples parameter ---
	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		bind_data->sample_indices =
		    ResolveSampleIndices(samples_it->second, bind_data->raw_sample_ct,
		                         bind_data->has_sample_info ? &bind_data->sample_info : nullptr, "read_pgen");

		bind_data->has_sample_subset = true;
		bind_data->subset_sample_ct = static_cast<uint32_t>(bind_data->sample_indices.size());
	}

	// --- Process variants parameter (single-file only — guarded above for multi-file) ---
	auto variants_it = input.named_parameters.find("variants");
	if (variants_it != input.named_parameters.end()) {
		auto resolved = ResolveVariantsParameter(variants_it->second, bind_data->Primary().variants,
		                                         bind_data->Primary().raw_variant_ct, "read_pgen");
		bind_data->has_variant_filter = true;
		bind_data->Primary().has_effective_variant_list = true;
		bind_data->Primary().effective_variant_indices = std::move(resolved);
	}

	// --- Cumulative effective-variant offsets across sources (global scan positions) ---
	bind_data->variant_offsets.push_back(0);
	{
		uint32_t cum = 0;
		for (auto &s : bind_data->sources) {
			cum += s.EffectiveVariantCt();
			bind_data->variant_offsets.push_back(cum);
		}
	}

	// --- Parse count filters (af_range, ac_range) ---
	{
		auto af_it = input.named_parameters.find("af_range");
		if (af_it != input.named_parameters.end()) {
			bind_data->count_filter.af_filter = ParseRangeFilter(af_it->second, "af_range", 0.0, 1.0, "read_pgen");
		}
		uint32_t pgen_output_sc = bind_data->has_sample_subset ? bind_data->subset_sample_ct : bind_data->sample_ct;
		auto ac_it = input.named_parameters.find("ac_range");
		if (ac_it != input.named_parameters.end()) {
			bind_data->count_filter.ac_filter =
			    ParseRangeFilter(ac_it->second, "ac_range", 0.0, static_cast<double>(2 * pgen_output_sc), "read_pgen");
		}

		if (bind_data->count_filter.HasFilter() && bind_data->has_sample_subset) {
			bind_data->count_filter_subset =
			    make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, bind_data->sample_indices));
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
			    "read_pgen: specify only one of include_genotypes or genotype_range (genotype_range is the numeric "
			    "alias of include_genotypes)");
		}
		if ((has_ig || has_gr) && bind_data->include_dosages) {
			throw InvalidInputException("read_pgen: %s is incompatible with dosages := true",
			                            has_ig ? "include_genotypes" : "genotype_range");
		}
		if (has_ig) {
			ParseIncludeGenotypes(ig_it->second, bind_data->genotype_filter, "read_pgen");
		} else if (has_gr) {
			bool inc_missing = false;
			RangeFilter range = ParseRangeFilter(gr_it->second, "genotype_range", 0.0, 2.0, "read_pgen", &inc_missing);
			bind_data->genotype_filter.SetFromRange(range, inc_missing);
		}
	}

	// --- Resolve genotype output mode ---
	uint32_t output_sample_ct = bind_data->has_sample_subset ? bind_data->subset_sample_ct : bind_data->sample_ct;

	string genotypes_str = "auto";
	auto genotypes_it = input.named_parameters.find("genotypes");
	if (genotypes_it != input.named_parameters.end()) {
		genotypes_str = genotypes_it->second.GetValue<string>();
	}
	bind_data->genotype_mode = ResolveGenotypeMode(genotypes_str, output_sample_ct, "read_pgen");

	// Validate incompatible combinations
	if (IsAggregateGenotypeMode(bind_data->genotype_mode)) {
		if (bind_data->include_phased) {
			throw InvalidInputException("read_pgen: genotypes := '%s' is incompatible with phased := true",
			                            bind_data->genotype_mode == GenotypeMode::COUNTS ? "counts" : "stats");
		}
		if (bind_data->include_dosages) {
			throw InvalidInputException("read_pgen: genotypes := '%s' is incompatible with dosages := true",
			                            bind_data->genotype_mode == GenotypeMode::COUNTS ? "counts" : "stats");
		}
	}

	// Build SampleSubset for COUNTS/STATS mode if needed
	if (IsAggregateGenotypeMode(bind_data->genotype_mode) && bind_data->has_sample_subset &&
	    !bind_data->count_filter_subset) {
		bind_data->count_filter_subset =
		    make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, bind_data->sample_indices));
	}

	// --- Register output columns ---
	if (bind_data->genotype_mode == GenotypeMode::COLUMNS) {
		// Columns mode: one scalar TINYINT column per output sample
		if (!bind_data->has_sample_info) {
			throw InvalidInputException("read_pgen: genotypes := 'columns' requires a .psam/.fam file for sample IDs "
			                            "(no companion file found)");
		}
		// Sort sample_indices for consistent pgenlib output order
		if (bind_data->has_sample_subset) {
			std::sort(bind_data->sample_indices.begin(), bind_data->sample_indices.end());
		}

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
		// STRUCT mode: single genotypes column with named fields per sample
		if (!bind_data->has_sample_info) {
			throw InvalidInputException("read_pgen: genotypes := 'struct' requires a .psam/.fam file for sample IDs "
			                            "(no companion file found)");
		}
		if (bind_data->has_sample_subset) {
			std::sort(bind_data->sample_indices.begin(), bind_data->sample_indices.end());
		}

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
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
		                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::STRUCT(std::move(struct_children))};
	} else if (bind_data->genotype_mode == GenotypeMode::COUNTS) {
		names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
		                LogicalType::VARCHAR, LogicalType::VARCHAR, MakeGenotypeCountsType()};
	} else if (bind_data->genotype_mode == GenotypeMode::STATS) {
		names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
		                LogicalType::VARCHAR, LogicalType::VARCHAR, MakeGenotypeStatsType()};
	} else {
		LogicalType elem_type = bind_data->include_phased    ? LogicalType::ARRAY(LogicalType::TINYINT, 2)
		                        : bind_data->include_dosages ? LogicalType::DOUBLE
		                                                     : LogicalType::TINYINT;
		LogicalType geno_type = bind_data->genotype_mode == GenotypeMode::ARRAY
		                            ? LogicalType::ARRAY(elem_type, output_sample_ct)
		                            : LogicalType::LIST(elem_type);

		names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
		                LogicalType::VARCHAR, LogicalType::VARCHAR, geno_type};
	}

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

//! Multi-file sub-batch size: no scan batch spans a file boundary, so a thread reopens
//! its reader at most once per claimed batch of this many variants.
static constexpr uint32_t PGEN_MULTI_SUBBATCH = 1000;

static unique_ptr<GlobalTableFunctionState> PgenInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PgenBindData>();
	auto state = make_uniq<PgenGlobalState>();

	state->total_variants = bind_data.EffectiveVariantCt();
	state->column_ids = input.column_ids;
	state->max_threads_config = GetPlinkingMaxThreads(context);

	// Check if genotypes column(s) are in the projection
	state->need_genotypes = false;
	bool need_aggregate_pgen = false;
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
		} else if (col_id == PgenBindData::GENOTYPES_COL_IDX) {
			if (IsAggregateGenotypeMode(bind_data.genotype_mode)) {
				// COUNTS/STATS use PgrGetCounts, not PgrGet
				need_aggregate_pgen = true;
			} else {
				state->need_genotypes = true;
			}
			break;
		}
	}

	state->need_pgen_reader = state->need_genotypes || bind_data.count_filter.HasFilter() ||
	                          bind_data.genotype_filter.active || need_aggregate_pgen;

	// Multi-file: precompute scan batches, each bounded to a single source and
	// sub-batched so no batch spans a file boundary. Threads claim batches atomically
	// (next_variant_idx as a batch index); a thread reopens its reader at most once per
	// claimed batch. Empty for single-file reads (they keep the fetch_add claim loop).
	if (bind_data.sources.size() > 1) {
		for (uint32_t s = 0; s < bind_data.sources.size(); s++) {
			uint32_t n = bind_data.sources[s].EffectiveVariantCt();
			for (uint32_t start = 0; start < n; start += PGEN_MULTI_SUBBATCH) {
				state->batches.push_back(ScanBatch {s, start, std::min(start + PGEN_MULTI_SUBBATCH, n)});
			}
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PgenInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                         GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PgenBindData>();
	auto &gstate = global_state->Cast<PgenGlobalState>();
	auto state = make_uniq<PgenLocalState>();

	if (!gstate.need_pgen_reader) {
		// No genotype columns or count filter needed — skip pgenlib initialization entirely
		return std::move(state);
	}

	// --- Allocate one-time per-thread buffers (sized off shared raw_sample_ct/sample_ct) ---
	// These are reused across all sources (identical samples by contract); only the
	// pgfi/pgr reader swaps at a file boundary (via OpenSourceReader).
	// Allocate genovec buffer (2 bits per sample, vector-aligned for SIMD safety)
	uint32_t effective_sample_ct = bind_data.has_sample_subset ? bind_data.raw_sample_ct : bind_data.sample_ct;
	uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(effective_sample_ct);
	state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));
	std::memset(state->genovec_buf.ptr, 0, genovec_word_ct * sizeof(uintptr_t));

	// Allocate expanded byte buffer for genotype conversion
	state->genotype_bytes.resize(effective_sample_ct);

	// Allocate phase buffers if phased output requested
	if (bind_data.include_phased) {
		uint32_t output_sample_ct = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.sample_ct;
		uintptr_t phase_wc = plink2::BitCtToAlignedWordCt(output_sample_ct);
		state->phasepresent_buf.Allocate(phase_wc * sizeof(uintptr_t));
		state->phaseinfo_buf.Allocate(phase_wc * sizeof(uintptr_t));
		state->phased_pairs.resize(static_cast<size_t>(output_sample_ct) * 2);
	}

	// Allocate dosage buffers if dosage output requested
	if (bind_data.include_dosages) {
		uint32_t output_sample_ct = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.sample_ct;
		uintptr_t dosage_present_wc = plink2::BitCtToAlignedWordCt(effective_sample_ct);
		state->dosage_present_buf.Allocate(dosage_present_wc * sizeof(uintptr_t));
		std::memset(state->dosage_present_buf.ptr, 0, dosage_present_wc * sizeof(uintptr_t));

		state->dosage_main_buf.Allocate(effective_sample_ct * sizeof(uint16_t));
		std::memset(state->dosage_main_buf.ptr, 0, effective_sample_ct * sizeof(uint16_t));

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

		// Build cumulative popcounts for PgrSampleSubsetIndex
		state->cumulative_popcounts_buf.Allocate(include_word_ct * sizeof(uint32_t));
		auto *cumulative_popcounts = state->cumulative_popcounts_buf.As<uint32_t>();
		plink2::FillCumulativePopcounts(sample_include, include_word_ct, cumulative_popcounts);
	}

	// Open the first source (single-file: the only one). Multi-file scans reopen on a
	// source boundary via OpenSourceReader. This binds pssi to the (new) pgr.
	OpenSourceReader(*state, bind_data, 0);

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static constexpr uint32_t PGEN_BATCH_SIZE = 128;

static void PgenScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PgenBindData>();
	auto &gstate = data_p.global_state->Cast<PgenGlobalState>();
	auto &lstate = data_p.local_state->Cast<PgenLocalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t total_variants = gstate.total_variants;

	// Effective sample count for output lists
	uint32_t output_sample_ct = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.sample_ct;

	idx_t rows_emitted = 0;

	// Emit one variant's row at output position `rows_emitted` (captured), reading
	// genotypes from `source` (the caller must have this source's reader open when
	// gstate.need_pgen_reader). Returns false if the variant was skipped by a count/
	// genotype pre-decompression filter (no row emitted). Shared by the single-file and
	// multi-file claim loops below.
	auto emit_variant_row = [&](const PgenSource &source, uint32_t local_pos) -> bool {
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
			uint32_t cf_sc = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.sample_ct;
			plink2::PglErr cf_err =
			    plink2::PgrGetCounts(cf_si, cf_iv, lstate.pssi, cf_sc, vidx, &lstate.pgr, genocounts);
			if (cf_err != plink2::kPglRetSuccess) {
				throw IOException("read_pgen: PgrGetCounts failed for variant %u", vidx);
			}
			auto pf = CheckPreDecompFilters(bind_data.count_filter, bind_data.genotype_filter, genocounts, cf_sc);
			if (pf.skip) {
				return false;
			}
			geno_range_all_pass = pf.all_pass;
		}

		// Read genotype data if needed (before filling columns, since
		// we need it for the genotypes column)
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
					throw IOException("read_pgen: PgrGetD failed for variant %u", vidx);
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
					throw IOException("read_pgen: PgrGetP failed for variant %u", vidx);
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
					throw IOException("read_pgen: PgrGet failed for variant %u", vidx);
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
			case 0: { // CHROM
				auto val = source.variants.GetChrom(vidx);
				FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				break;
			}
			case 1: { // POS
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = source.variants.GetPos(vidx);
				break;
			}
			case 2: { // ID
				auto val = source.variants.GetId(vidx);
				if (val.empty()) {
					FlatVector::SetNull(vec, rows_emitted, true);
				} else {
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				}
				break;
			}
			case 3: { // REF
				auto val = source.variants.GetRef(vidx);
				FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				break;
			}
			case 4: { // ALT
				auto val = source.variants.GetAlt(vidx);
				if (val.empty() || val == ".") {
					FlatVector::SetNull(vec, rows_emitted, true);
				} else {
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
				}
				break;
			}
			case 5: { // genotypes
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
					uint32_t agg_sc = bind_data.has_sample_subset ? bind_data.subset_sample_ct : bind_data.sample_ct;

					plink2::PglErr agg_err =
					    plink2::PgrGetCounts(agg_si, agg_iv, lstate.pssi, agg_sc, vidx, &lstate.pgr, genocounts);
					if (agg_err != plink2::kPglRetSuccess) {
						throw IOException("read_pgen: PgrGetCounts failed for variant %u", vidx);
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
					// Phased output: ARRAY(ARRAY(TINYINT,2), N) or LIST(ARRAY(TINYINT,2))
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
							if (a1 == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
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
						// LIST(ARRAY(TINYINT, 2))
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
							if (a1 == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
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
				} else {
					// Unphased output: ARRAY(TINYINT, N) or LIST(TINYINT)
					if (bind_data.genotype_mode == GenotypeMode::ARRAY) {
						auto array_size = static_cast<idx_t>(output_sample_ct);
						auto &child = ArrayVector::GetEntry(vec);
						auto *child_data = FlatVector::GetData<int8_t>(child);
						auto &child_validity = FlatVector::Validity(child);

						idx_t base = rows_emitted * array_size;
						for (idx_t s = 0; s < array_size; s++) {
							int8_t geno = lstate.genotype_bytes[s];
							if (geno == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
							                   !bind_data.genotype_filter.AllowsCall(static_cast<double>(geno)))) {
								child_validity.SetInvalid(base + s);
								child_data[base + s] = 0;
							} else {
								child_data[base + s] = geno;
							}
						}
					} else {
						// LIST path
						auto list_offset = ListVector::GetListSize(vec);
						ListVector::Reserve(vec, list_offset + output_sample_ct);
						auto &child = ListVector::GetEntry(vec);
						auto *child_data = FlatVector::GetData<int8_t>(child);
						auto &child_validity = FlatVector::Validity(child);
						for (idx_t s = 0; s < output_sample_ct; s++) {
							int8_t geno = lstate.genotype_bytes[s];
							if (geno == -9 || (bind_data.genotype_filter.active && !geno_range_all_pass &&
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
		const PgenSource &source = bind_data.Primary();
		while (rows_emitted < STANDARD_VECTOR_SIZE) {
			// Claim a batch of variants, capped to remaining output capacity
			uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
			uint32_t claim_size = std::min(PGEN_BATCH_SIZE, remaining_capacity);
			uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
			if (batch_start >= total_variants) {
				break;
			}
			uint32_t batch_end = std::min(batch_start + claim_size, total_variants);
			for (uint32_t ev = batch_start; ev < batch_end; ev++) {
				if (emit_variant_row(source, ev)) {
					rows_emitted++;
				}
			}
		}
	} else {
		// Multi-file: claim precomputed batches (each bounded to one source), reopening the
		// per-thread reader only at a source boundary (OpenSourceReader is a no-op otherwise).
		while (rows_emitted < STANDARD_VECTOR_SIZE) {
			if (lstate.mf_need_claim) {
				lstate.mf_batch = gstate.next_variant_idx.fetch_add(1);
				if (lstate.mf_batch >= gstate.batches.size()) {
					break;
				}
				lstate.mf_local = gstate.batches[lstate.mf_batch].local_start;
				lstate.mf_need_claim = false;
			}
			const ScanBatch &batch = gstate.batches[lstate.mf_batch];
			const PgenSource &source = bind_data.sources[batch.source_idx];
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
// Registration
// ---------------------------------------------------------------------------

void RegisterPgenReader(ExtensionLoader &loader) {
	// Two overloads: a single .pgen path (VARCHAR) and a list of .pgen paths
	// (LIST(VARCHAR), row-concatenated across variant-sharded pgens). Both dispatch to
	// the same callbacks; bind detects the list via ResolvePathList.
	auto add_named_params = [](TableFunction &fn) {
		fn.projection_pushdown = true;
		fn.named_parameters["pvar"] = LogicalType::VARCHAR;
		fn.named_parameters["psam"] = LogicalType::VARCHAR;
		fn.named_parameters["dosages"] = LogicalType::BOOLEAN;
		fn.named_parameters["phased"] = LogicalType::BOOLEAN;
		// Accept ANY for samples — type dispatch (LIST(INTEGER) vs LIST(VARCHAR))
		// is handled in PgenBind based on the actual value type.
		fn.named_parameters["samples"] = LogicalType::ANY;
		fn.named_parameters["genotypes"] = LogicalType::VARCHAR;
		fn.named_parameters["orient"] = LogicalType::VARCHAR;
		fn.named_parameters["af_range"] = LogicalType::ANY;
		fn.named_parameters["ac_range"] = LogicalType::ANY;
		fn.named_parameters["genotype_range"] = LogicalType::ANY;
		fn.named_parameters["include_genotypes"] = LogicalType::LIST(LogicalType::VARCHAR);
		fn.named_parameters["variants"] = LogicalType::ANY;
	};

	TableFunctionSet set("read_pgen");
	TableFunction one("read_pgen", {LogicalType::VARCHAR}, PgenScan, PgenBind, PgenInitGlobal, PgenInitLocal);
	add_named_params(one);
	set.AddFunction(one);
	TableFunction many("read_pgen", {LogicalType::LIST(LogicalType::VARCHAR)}, PgenScan, PgenBind, PgenInitGlobal,
	                   PgenInitLocal);
	add_named_params(many);
	set.AddFunction(many);
	loader.RegisterFunction(set);
}

} // namespace duckdb
