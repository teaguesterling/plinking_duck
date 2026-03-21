#include "plink_pca.hpp"
#include "plink_common.hpp"

#include "duckdb/parallel/task_scheduler.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <mutex>
#include <random>

namespace duckdb {

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

static constexpr uint32_t PCA_VARIANT_BLOCK_SIZE = 240;

// ---------------------------------------------------------------------------
// Output modes
// ---------------------------------------------------------------------------

enum class PcaMode : uint8_t { SAMPLES, PCS, BOTH };

static PcaMode ParsePcaMode(const string &mode_str) {
	auto lower = StringUtil::Lower(mode_str);
	if (lower == "samples") {
		return PcaMode::SAMPLES;
	} else if (lower == "pcs") {
		return PcaMode::PCS;
	} else if (lower == "both") {
		return PcaMode::BOTH;
	}
	throw InvalidInputException("plink_pca: invalid mode '%s' (expected 'samples', 'pcs', or 'both')", mode_str);
}

// ---------------------------------------------------------------------------
// Column indices — samples mode
// ---------------------------------------------------------------------------

static constexpr idx_t SCOL_FID = 0;
static constexpr idx_t SCOL_IID = 1;
static constexpr idx_t SCOL_PC_START = 2; // PC1..PCk start at column 2

// ---------------------------------------------------------------------------
// Column indices — pcs mode
// ---------------------------------------------------------------------------

static constexpr idx_t PCOL_PC = 0;
static constexpr idx_t PCOL_EIGENVALUE = 1;
static constexpr idx_t PCOL_VARIANCE_PROPORTION = 2;
static constexpr idx_t PCOL_CUMULATIVE_VARIANCE = 3;

// ---------------------------------------------------------------------------
// Effective variant (pgen index + normalization)
// ---------------------------------------------------------------------------

struct EffectiveVariant {
	uint32_t pgen_idx;
	VariantNorm norm;
};

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkPcaBindData : public TableFunctionData {
	string pgen_path;
	string pvar_path;
	string psam_path;

	VariantMetadataIndex variants;
	SampleInfo sample_info;

	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Sample subsetting
	bool has_sample_subset = false;
	unique_ptr<SampleSubset> sample_subset;
	uint32_t effective_sample_ct = 0;

	// Region filtering
	VariantRange variant_range;

	// Effective variants with pre-computed normalization
	vector<EffectiveVariant> effective_variants;
	uint32_t effective_variant_ct = 0;

	// Output
	PcaMode mode = PcaMode::SAMPLES;
	uint32_t n_pcs = 10;

	// Derived algorithm constants
	uint32_t pc_ct_x2 = 0;  // 2 * n_pcs
	uint32_t qq_col_ct = 0; // (n_pcs + 1) * pc_ct_x2

	// Sample output order (maps emit index → original sample index)
	vector<uint32_t> sample_output_order;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkPcaGlobalState : public GlobalTableFunctionState {
	// Algorithm dimensions
	uint32_t N = 0; // effective sample count
	uint32_t M = 0; // effective variant count
	uint32_t n_pcs = 0;
	uint32_t pc_ct_x2 = 0;
	uint32_t qq_col_ct = 0;

	// Shared algorithm matrices (row-major)
	vector<double> G1; // N x pc_ct_x2
	vector<double> QQ; // M x qq_col_ct

	// Per-thread partial accumulation buffers
	// thread_partials[tid] is N x qq_col_ct doubles
	uint32_t thread_count = 0;
	vector<vector<double>> thread_partials;

	// Pass coordination — generation-based barrier
	//
	// All threads participate in every pass. Each thread increments
	// pass_active_threads on entry and decrements on finish. The last
	// thread to finish (decrement returns 1) does the merge/SVD, resets
	// next_block_idx, and increments pass_generation to release waiters.
	// Non-last threads return empty chunks and check pass_generation on
	// re-entry — if it has advanced, they join the new pass.
	uint32_t total_passes = 0;                     // n_pcs + 2
	std::atomic<uint32_t> pass_generation {0};     // incremented after each pass completes
	std::atomic<uint32_t> pass_active_threads {0}; // threads currently in a pass
	std::atomic<uint32_t> next_block_idx {0};      // variant block claiming
	std::atomic<bool> algorithm_done {false};

	// Thread ID assignment
	std::atomic<uint32_t> next_thread_id {0};

	// Results
	vector<double> eigenvectors; // N x n_pcs (row-major)
	vector<double> eigenvalues;  // n_pcs

	// Emission
	std::atomic<uint32_t> next_emit_idx {0};
	vector<column_t> column_ids;

	// Threading
	uint32_t db_thread_count = 1;

	idx_t MaxThreads() const override {
		if (M < PCA_VARIANT_BLOCK_SIZE) {
			return 1;
		}
		return std::min<idx_t>(M / PCA_VARIANT_BLOCK_SIZE + 1, db_thread_count);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkPcaLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	AlignedBuffer genovec_buf;
	vector<int8_t> geno_bytes;
	vector<double> norm_geno;

	uint32_t thread_id = 0;
	uint32_t last_generation_seen = UINT32_MAX;
	bool initialized = false;

	~PlinkPcaLocalState() {
		if (initialized) {
			plink2::PglErr reterr = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr, &reterr);
			plink2::CleanupPgfi(&pgfi, &reterr);
		}
	}
};

// ---------------------------------------------------------------------------
// Bind
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PlinkPcaBind(ClientContext &context, TableFunctionBindInput &input,
                                             vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkPcaBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters (first pass) ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "mode") {
			bind_data->mode = ParsePcaMode(kv.second.GetValue<string>());
		} else if (kv.first == "n_pcs") {
			int32_t val = kv.second.GetValue<int32_t>();
			if (val < 1) {
				throw InvalidInputException("plink_pca: n_pcs must be >= 1 (got %d)", val);
			}
			bind_data->n_pcs = static_cast<uint32_t>(val);
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled below
		}
	}

	// Compute algorithm constants
	bind_data->pc_ct_x2 = 2 * bind_data->n_pcs;
	bind_data->qq_col_ct = (bind_data->n_pcs + 1) * bind_data->pc_ct_x2;

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_pca: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		if (bind_data->psam_path.empty()) {
			throw InvalidInputException("plink_pca: cannot find .psam or .fam companion for '%s' "
			                            "(use psam := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	// --- Initialize pgenlib (temporary, for header + allele freq counting) ---
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
		throw IOException("plink_pca: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;

	AlignedBuffer pgfi_alloc;
	if (pgfi_alloc_cacheline_ct > 0) {
		pgfi_alloc.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	}

	uint32_t max_vrec_width = 0;
	uintptr_t pgr_alloc_cacheline_ct = 0;

	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, pgfi.raw_variant_ct, &max_vrec_width, &pgfi,
	                             pgfi_alloc.As<unsigned char>(), &pgr_alloc_cacheline_ct, errstr_buf);
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("plink_pca: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadataIndex(context, bind_data->pvar_path, "plink_pca");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw InvalidInputException("plink_pca: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info ---
	bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);

	if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw InvalidInputException("plink_pca: sample count mismatch: .pgen has %u samples, "
		                            ".psam/.fam '%s' has %llu samples",
		                            bind_data->raw_sample_ct, bind_data->psam_path,
		                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
	}

	// --- Process samples parameter ---
	bind_data->effective_sample_ct = bind_data->raw_sample_ct;

	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		auto indices =
		    ResolveSampleIndices(samples_it->second, bind_data->raw_sample_ct, &bind_data->sample_info, "plink_pca");

		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;

		auto sorted_indices = indices;
		std::sort(sorted_indices.begin(), sorted_indices.end());
		bind_data->sample_output_order = std::move(sorted_indices);
	} else {
		bind_data->sample_output_order.resize(bind_data->raw_sample_ct);
		for (uint32_t i = 0; i < bind_data->raw_sample_ct; i++) {
			bind_data->sample_output_order[i] = i;
		}
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range = ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_pca");
	}

	// --- Compute per-variant allele frequencies and build effective variant list ---
	// Need a full PgenReader for PgrGetCounts
	plink2::PgenReader pgr_temp;
	plink2::PreinitPgr(&pgr_temp);
	AlignedBuffer pgr_alloc_temp;
	if (pgr_alloc_cacheline_ct > 0) {
		pgr_alloc_temp.Allocate(pgr_alloc_cacheline_ct * plink2::kCacheline);
	}

	err = plink2::PgrInit(bind_data->pgen_path.c_str(), max_vrec_width, &pgfi, &pgr_temp,
	                      pgr_alloc_temp.As<unsigned char>());
	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&pgr_temp, &cleanup_err);
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("plink_pca: PgrInit failed for '%s'", bind_data->pgen_path);
	}

	// Set up sample subsetting on the temporary reader
	plink2::PgrSampleSubsetIndex pssi_temp;
	const uintptr_t *sample_include_temp = nullptr;
	const uintptr_t *interleaved_vec_temp = nullptr;
	uint32_t count_sample_ct = bind_data->effective_sample_ct;

	if (bind_data->has_sample_subset && bind_data->sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data->sample_subset->CumulativePopcounts(), &pgr_temp, &pssi_temp);
		sample_include_temp = bind_data->sample_subset->SampleInclude();
		interleaved_vec_temp = bind_data->sample_subset->InterleavedVec();
	} else {
		plink2::PgrClearSampleSubsetIndex(&pgr_temp, &pssi_temp);
	}

	uint32_t range_start = bind_data->variant_range.has_filter ? bind_data->variant_range.start_idx : 0;
	uint32_t range_end =
	    bind_data->variant_range.has_filter ? bind_data->variant_range.end_idx : bind_data->raw_variant_ct;

	for (uint32_t vidx = range_start; vidx < range_end; vidx++) {
		STD_ARRAY_DECL(uint32_t, 4, genocounts);
		err = plink2::PgrGetCounts(sample_include_temp, interleaved_vec_temp, pssi_temp, count_sample_ct, vidx,
		                           &pgr_temp, genocounts);
		if (err != plink2::kPglRetSuccess) {
			plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr_temp, &cleanup_err);
			plink2::CleanupPgfi(&pgfi, &cleanup_err);
			throw IOException("plink_pca: PgrGetCounts failed for variant %u", vidx);
		}

		uint32_t obs = genocounts[0] + genocounts[1] + genocounts[2];
		if (obs == 0) {
			continue; // all missing
		}

		double alt_freq = (static_cast<double>(genocounts[1]) + 2.0 * static_cast<double>(genocounts[2])) /
		                  (2.0 * static_cast<double>(obs));
		VariantNorm norm = ComputeVariantNorm(alt_freq);
		if (norm.skip) {
			continue; // monomorphic
		}

		bind_data->effective_variants.push_back({vidx, norm});
	}

	// Clean up temporary reader
	{
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&pgr_temp, &cleanup_err);
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
	}

	bind_data->effective_variant_ct = static_cast<uint32_t>(bind_data->effective_variants.size());

	// --- Validate dimensions ---
	if (bind_data->effective_sample_ct < 2) {
		throw InvalidInputException("plink_pca: need at least 2 samples (got %u)", bind_data->effective_sample_ct);
	}

	if (bind_data->n_pcs >= bind_data->effective_sample_ct) {
		throw InvalidInputException("plink_pca: n_pcs (%u) must be less than sample count (%u)", bind_data->n_pcs,
		                            bind_data->effective_sample_ct);
	}

	if (bind_data->effective_variant_ct <= bind_data->qq_col_ct) {
		throw InvalidInputException(
		    "plink_pca: too few variants (%u) for %u PCs with approx mode (need > %u non-monomorphic variants)",
		    bind_data->effective_variant_ct, bind_data->n_pcs, bind_data->qq_col_ct);
	}

	if (bind_data->effective_sample_ct <= bind_data->qq_col_ct) {
		throw InvalidInputException("plink_pca: too few samples (%u) for %u PCs with approx mode "
		                            "(need > %u samples; try fewer PCs or more samples)",
		                            bind_data->effective_sample_ct, bind_data->n_pcs, bind_data->qq_col_ct);
	}

	// --- Memory guard ---
	uint64_t qq_elements = static_cast<uint64_t>(bind_data->effective_variant_ct) * bind_data->qq_col_ct;
	Value max_elements_val;
	uint64_t max_elements = 16ULL * 1024 * 1024 * 1024; // default
	if (context.TryGetCurrentSetting("plinking_max_matrix_elements", max_elements_val)) {
		auto val = max_elements_val.GetValue<int64_t>();
		max_elements = val > 0 ? static_cast<uint64_t>(val) : 0;
	}
	if (qq_elements > max_elements) {
		throw InvalidInputException(
		    "plink_pca: QQ matrix would require %llu elements (%llu MB), exceeding plinking_max_matrix_elements "
		    "(%lld). "
		    "Reduce n_pcs or variant count, or increase the limit with SET plinking_max_matrix_elements = <value>.",
		    static_cast<unsigned long long>(qq_elements),
		    static_cast<unsigned long long>(qq_elements * 8 / (1024 * 1024)), static_cast<long long>(max_elements));
	}

	// --- Register output schema ---
	if (bind_data->mode == PcaMode::PCS) {
		names = {"PC", "EIGENVALUE", "VARIANCE_PROPORTION", "CUMULATIVE_VARIANCE"};
		return_types = {LogicalType::INTEGER, LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE};
	} else if (bind_data->mode == PcaMode::SAMPLES) {
		names = {"FID", "IID"};
		return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR};
		for (uint32_t i = 0; i < bind_data->n_pcs; i++) {
			names.push_back("PC" + std::to_string(i + 1));
			return_types.push_back(LogicalType::DOUBLE);
		}
	} else {
		// mode = 'both': nested struct
		// EIGENVEC: LIST(STRUCT(FID, IID, PC1..PCk))
		// EIGENVAL: LIST(DOUBLE)
		child_list_t<LogicalType> eigenvec_fields;
		eigenvec_fields.push_back({"FID", LogicalType::VARCHAR});
		eigenvec_fields.push_back({"IID", LogicalType::VARCHAR});
		for (uint32_t i = 0; i < bind_data->n_pcs; i++) {
			eigenvec_fields.push_back({"PC" + std::to_string(i + 1), LogicalType::DOUBLE});
		}
		auto eigenvec_struct = LogicalType::STRUCT(std::move(eigenvec_fields));
		auto eigenvec_list = LogicalType::LIST(std::move(eigenvec_struct));

		names = {"EIGENVEC", "EIGENVAL"};
		return_types = {std::move(eigenvec_list), LogicalType::LIST(LogicalType::DOUBLE)};
	}

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkPcaInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkPcaBindData>();
	auto state = make_uniq<PlinkPcaGlobalState>();

	state->N = bind_data.effective_sample_ct;
	state->M = bind_data.effective_variant_ct;
	state->n_pcs = bind_data.n_pcs;
	state->pc_ct_x2 = bind_data.pc_ct_x2;
	state->qq_col_ct = bind_data.qq_col_ct;
	state->total_passes = bind_data.n_pcs + 2;

	state->column_ids = input.column_ids;
	state->db_thread_count = static_cast<uint32_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	state->thread_count = static_cast<uint32_t>(state->MaxThreads());

	// Initialize G1 with Gaussian random noise (fixed seed)
	state->G1.resize(static_cast<size_t>(state->N) * state->pc_ct_x2);
	std::mt19937_64 rng(12345);
	std::normal_distribution<double> dist(0.0, 1.0);
	for (auto &val : state->G1) {
		val = dist(rng);
	}

	// Allocate QQ (M x qq_col_ct)
	state->QQ.resize(static_cast<size_t>(state->M) * state->qq_col_ct, 0.0);

	// Allocate per-thread partial buffers
	state->thread_partials.resize(state->thread_count);
	for (auto &partial : state->thread_partials) {
		partial.resize(static_cast<size_t>(state->N) * state->qq_col_ct, 0.0);
	}

	// Allocate results
	state->eigenvectors.resize(static_cast<size_t>(state->N) * state->n_pcs, 0.0);
	state->eigenvalues.resize(state->n_pcs, 0.0);

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkPcaInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                             GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkPcaBindData>();
	auto &gstate = global_state->Cast<PlinkPcaGlobalState>();
	auto state = make_uniq<PlinkPcaLocalState>();

	if (bind_data.effective_variants.empty()) {
		return std::move(state);
	}

	// Assign thread ID
	state->thread_id = gstate.next_thread_id.fetch_add(1);

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
		throw IOException("plink_pca: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_pca: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_pca: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	// Allocate genotype decode buffer
	uint32_t raw_sample_ct = bind_data.raw_sample_ct;
	uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(raw_sample_ct);
	state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));
	std::memset(state->genovec_buf.ptr, 0, genovec_word_ct * sizeof(uintptr_t));

	// Allocate genotype processing buffers
	state->geno_bytes.resize(bind_data.effective_sample_ct);
	state->norm_geno.resize(bind_data.effective_sample_ct);

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Accumulation helpers
// ---------------------------------------------------------------------------

// Step A: qq_row = norm_geno (1 x N) × G1 (N x pc_ct_x2) → 1 x pc_ct_x2
// Written into QQ at row eff_idx, columns [col_offset, col_offset + pc_ct_x2)
static void AccumulateStepA(PlinkPcaGlobalState &gs, const double *norm_geno, uint32_t eff_idx, uint32_t col_offset) {
	uint32_t N = gs.N;
	uint32_t pc_ct_x2 = gs.pc_ct_x2;
	double *qq_row = &gs.QQ[static_cast<size_t>(eff_idx) * gs.qq_col_ct + col_offset];
	const double *g1 = gs.G1.data();

	for (uint32_t c = 0; c < pc_ct_x2; c++) {
		double sum = 0.0;
		for (uint32_t s = 0; s < N; s++) {
			sum += norm_geno[s] * g1[static_cast<size_t>(s) * pc_ct_x2 + c];
		}
		qq_row[c] = sum;
	}
}

// Step B: g2_part += norm_genoᵀ (N x 1) × qq_row (1 x pc_ct_x2) → rank-1 update N x pc_ct_x2
// Accumulated into thread_partials[tid] at appropriate column offset
static void AccumulateStepB(PlinkPcaGlobalState &gs, uint32_t tid, const double *norm_geno, const double *qq_row,
                            uint32_t col_offset) {
	uint32_t N = gs.N;
	uint32_t pc_ct_x2 = gs.pc_ct_x2;
	double *partial = gs.thread_partials[tid].data();

	for (uint32_t s = 0; s < N; s++) {
		double g = norm_geno[s];
		for (uint32_t c = 0; c < pc_ct_x2; c++) {
			partial[static_cast<size_t>(s) * gs.qq_col_ct + col_offset + c] += g * qq_row[c];
		}
	}
}

// Phase 3: bb_part += norm_genoᵀ (N x 1) × qq_row (1 x qq_col_ct) → rank-1 update N x qq_col_ct
static void AccumulatePhase3(PlinkPcaGlobalState &gs, uint32_t tid, const double *norm_geno, uint32_t eff_idx) {
	uint32_t N = gs.N;
	uint32_t qq_col_ct = gs.qq_col_ct;
	const double *qq_row = &gs.QQ[static_cast<size_t>(eff_idx) * qq_col_ct];
	double *partial = gs.thread_partials[tid].data();

	for (uint32_t s = 0; s < N; s++) {
		double g = norm_geno[s];
		for (uint32_t c = 0; c < qq_col_ct; c++) {
			partial[static_cast<size_t>(s) * qq_col_ct + c] += g * qq_row[c];
		}
	}
}

// ---------------------------------------------------------------------------
// SVD helpers (single-threaded, Eigen3)
// ---------------------------------------------------------------------------

// Phase 2 SVD: QQ (M x qq_col_ct) → overwrite with left singular vectors
static void RunKrylovSVD(PlinkPcaGlobalState &gs) {
	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> qq_map(gs.QQ.data(), gs.M,
	                                                                                          gs.qq_col_ct);

	Eigen::BDCSVD<Eigen::MatrixXd> svd(qq_map, Eigen::ComputeThinU);

	// Overwrite QQ with left singular vectors (M x qq_col_ct)
	// Eigen BDCSVD returns column-major U; we need to write back row-major
	auto &U = svd.matrixU();
	for (uint32_t r = 0; r < gs.M; r++) {
		for (uint32_t c = 0; c < gs.qq_col_ct; c++) {
			gs.QQ[static_cast<size_t>(r) * gs.qq_col_ct + c] = U(r, c);
		}
	}
}

// Phase 3 SVD: BB (N x qq_col_ct) → eigenvectors and eigenvalues
static void RunFinalSVD(PlinkPcaGlobalState &gs, const vector<double> &bb) {
	Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> bb_map(bb.data(), gs.N,
	                                                                                                gs.qq_col_ct);

	Eigen::BDCSVD<Eigen::MatrixXd> svd(bb_map, Eigen::ComputeThinU);

	auto &U = svd.matrixU();
	auto &S = svd.singularValues();

	// Extract top n_pcs eigenvectors (columns of U) and eigenvalues (S²/M)
	for (uint32_t s = 0; s < gs.N; s++) {
		for (uint32_t pc = 0; pc < gs.n_pcs; pc++) {
			gs.eigenvectors[static_cast<size_t>(s) * gs.n_pcs + pc] = U(s, pc);
		}
	}

	double total_variance = 0.0;
	for (uint32_t pc = 0; pc < gs.n_pcs; pc++) {
		gs.eigenvalues[pc] = S(pc) * S(pc) / static_cast<double>(gs.M);
		total_variance += gs.eigenvalues[pc];
	}

	// Store total variance for proportion computation (store after eigenvalues)
	// We'll compute proportions at emission time using the eigenvalues directly
}

// ---------------------------------------------------------------------------
// Emission helpers
// ---------------------------------------------------------------------------

static void EmitSamplesMode(const PlinkPcaBindData &bind_data, PlinkPcaGlobalState &gs, DataChunk &output) {
	auto &column_ids = gs.column_ids;
	bool has_fid = !bind_data.sample_info.fids.empty();
	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t sidx = gs.next_emit_idx.fetch_add(1);
		if (sidx >= gs.N) {
			break;
		}

		uint32_t orig_idx = bind_data.sample_output_order[sidx];

		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];
			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];

			if (file_col == SCOL_FID) {
				if (has_fid) {
					FlatVector::GetData<string_t>(vec)[rows_emitted] =
					    StringVector::AddString(vec, bind_data.sample_info.fids[orig_idx]);
				} else {
					FlatVector::SetNull(vec, rows_emitted, true);
				}
			} else if (file_col == SCOL_IID) {
				FlatVector::GetData<string_t>(vec)[rows_emitted] =
				    StringVector::AddString(vec, bind_data.sample_info.iids[orig_idx]);
			} else if (file_col >= SCOL_PC_START && file_col < SCOL_PC_START + bind_data.n_pcs) {
				uint32_t pc = static_cast<uint32_t>(file_col - SCOL_PC_START);
				FlatVector::GetData<double>(vec)[rows_emitted] =
				    gs.eigenvectors[static_cast<size_t>(sidx) * gs.n_pcs + pc];
			}
		}
		rows_emitted++;
	}
	output.SetCardinality(rows_emitted);
}

static void EmitPcsMode(const PlinkPcaBindData &bind_data, PlinkPcaGlobalState &gs, DataChunk &output) {
	auto &column_ids = gs.column_ids;
	idx_t rows_emitted = 0;

	// Compute total variance for proportions
	double total_variance = 0.0;
	for (uint32_t i = 0; i < gs.n_pcs; i++) {
		total_variance += gs.eigenvalues[i];
	}

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t pc_idx = gs.next_emit_idx.fetch_add(1);
		if (pc_idx >= gs.n_pcs) {
			break;
		}

		double eigenvalue = gs.eigenvalues[pc_idx];
		double var_prop = (total_variance > 0.0) ? eigenvalue / total_variance : 0.0;

		// Compute cumulative variance
		double cum_var = 0.0;
		for (uint32_t i = 0; i <= pc_idx; i++) {
			cum_var += gs.eigenvalues[i];
		}
		cum_var = (total_variance > 0.0) ? cum_var / total_variance : 0.0;

		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];
			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];

			switch (file_col) {
			case PCOL_PC:
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(pc_idx + 1);
				break;
			case PCOL_EIGENVALUE:
				FlatVector::GetData<double>(vec)[rows_emitted] = eigenvalue;
				break;
			case PCOL_VARIANCE_PROPORTION:
				FlatVector::GetData<double>(vec)[rows_emitted] = var_prop;
				break;
			case PCOL_CUMULATIVE_VARIANCE:
				FlatVector::GetData<double>(vec)[rows_emitted] = cum_var;
				break;
			default:
				break;
			}
		}
		rows_emitted++;
	}
	output.SetCardinality(rows_emitted);
}

static void EmitBothMode(const PlinkPcaBindData &bind_data, PlinkPcaGlobalState &gs, DataChunk &output) {
	// Single row with nested structures
	uint32_t emit_idx = gs.next_emit_idx.fetch_add(1);
	if (emit_idx > 0) {
		output.SetCardinality(0);
		return;
	}

	bool has_fid = !bind_data.sample_info.fids.empty();

	// Build EIGENVEC list: LIST(STRUCT(FID, IID, PC1..PCk))
	vector<Value> eigenvec_entries;
	for (uint32_t sidx = 0; sidx < gs.N; sidx++) {
		uint32_t orig_idx = bind_data.sample_output_order[sidx];

		child_list_t<Value> fields;
		fields.push_back({"FID", has_fid ? Value(bind_data.sample_info.fids[orig_idx]) : Value()});
		fields.push_back({"IID", Value(bind_data.sample_info.iids[orig_idx])});
		for (uint32_t pc = 0; pc < gs.n_pcs; pc++) {
			fields.push_back({"PC" + std::to_string(pc + 1),
			                  Value::DOUBLE(gs.eigenvectors[static_cast<size_t>(sidx) * gs.n_pcs + pc])});
		}
		eigenvec_entries.push_back(Value::STRUCT(std::move(fields)));
	}

	// Build EIGENVAL list: LIST(DOUBLE)
	vector<Value> eigenval_entries;
	for (uint32_t pc = 0; pc < gs.n_pcs; pc++) {
		eigenval_entries.push_back(Value::DOUBLE(gs.eigenvalues[pc]));
	}

	auto &column_ids = gs.column_ids;
	for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
		auto file_col = column_ids[out_col];
		if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}

		auto &vec = output.data[out_col];

		if (file_col == 0) {
			vec.SetValue(0, Value::LIST(std::move(eigenvec_entries)));
		} else if (file_col == 1) {
			vec.SetValue(0, Value::LIST(std::move(eigenval_entries)));
		}
	}
	output.SetCardinality(1);
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void ScanVariantPass(const PlinkPcaBindData &bind_data, PlinkPcaGlobalState &gs, PlinkPcaLocalState &ls,
                            const uintptr_t *sample_include, uint32_t pass) {
	uint32_t sample_ct = bind_data.effective_sample_ct;
	uint32_t M = gs.M;
	bool is_phase3 = (pass == bind_data.n_pcs + 1);
	uint32_t col_offset = is_phase3 ? 0 : pass * gs.pc_ct_x2;

	while (true) {
		uint32_t block_start = gs.next_block_idx.fetch_add(PCA_VARIANT_BLOCK_SIZE);
		if (block_start >= M) {
			break;
		}
		uint32_t block_end = std::min(block_start + PCA_VARIANT_BLOCK_SIZE, M);

		for (uint32_t eff_idx = block_start; eff_idx < block_end; eff_idx++) {
			auto &ev = bind_data.effective_variants[eff_idx];

			plink2::PglErr err = plink2::PgrGet(sample_include, ls.pssi, sample_ct, ev.pgen_idx, &ls.pgr,
			                                    ls.genovec_buf.As<uintptr_t>());
			if (err != plink2::kPglRetSuccess) {
				throw IOException("plink_pca: PgrGet failed for variant %u", ev.pgen_idx);
			}

			plink2::GenoarrToBytesMinus9(ls.genovec_buf.As<uintptr_t>(), sample_ct, ls.geno_bytes.data());
			NormalizeGenotypes(ls.geno_bytes.data(), sample_ct, ev.norm, ls.norm_geno.data());

			if (is_phase3) {
				AccumulatePhase3(gs, ls.thread_id, ls.norm_geno.data(), eff_idx);
			} else {
				AccumulateStepA(gs, ls.norm_geno.data(), eff_idx, col_offset);

				if (pass < bind_data.n_pcs) {
					const double *qq_row = &gs.QQ[static_cast<size_t>(eff_idx) * gs.qq_col_ct + col_offset];
					AccumulateStepB(gs, ls.thread_id, ls.norm_geno.data(), qq_row, col_offset);
				}
			}
		}
	}
}

static void MergePass(const PlinkPcaBindData &bind_data, PlinkPcaGlobalState &gs, uint32_t pass) {
	uint32_t N = gs.N;
	bool is_phase3 = (pass == bind_data.n_pcs + 1);
	uint32_t col_offset = is_phase3 ? 0 : pass * gs.pc_ct_x2;

	if (is_phase3) {
		vector<double> bb(static_cast<size_t>(N) * gs.qq_col_ct, 0.0);
		for (uint32_t t = 0; t < gs.thread_count; t++) {
			auto &partial = gs.thread_partials[t];
			for (size_t i = 0; i < bb.size(); i++) {
				bb[i] += partial[i];
			}
		}
		RunFinalSVD(gs, bb);
		gs.algorithm_done.store(true, std::memory_order_release);
		return;
	}

	// Merge thread partials → G2, compute G1 = G2 / M
	std::fill(gs.G1.begin(), gs.G1.end(), 0.0);
	for (uint32_t t = 0; t < gs.thread_count; t++) {
		auto &partial = gs.thread_partials[t];
		for (uint32_t s = 0; s < N; s++) {
			for (uint32_t c = 0; c < gs.pc_ct_x2; c++) {
				gs.G1[static_cast<size_t>(s) * gs.pc_ct_x2 + c] +=
				    partial[static_cast<size_t>(s) * gs.qq_col_ct + col_offset + c];
			}
		}
	}

	double inv_M = 1.0 / static_cast<double>(gs.M);
	for (auto &val : gs.G1) {
		val *= inv_M;
	}

	if (pass == bind_data.n_pcs) {
		RunKrylovSVD(gs);
	}
}

static void PlinkPcaScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkPcaBindData>();
	auto &gs = data_p.global_state->Cast<PlinkPcaGlobalState>();
	auto &ls = data_p.local_state->Cast<PlinkPcaLocalState>();

	// --- Emit rows (algorithm complete) ---
	if (gs.algorithm_done.load(std::memory_order_acquire)) {
		switch (bind_data.mode) {
		case PcaMode::SAMPLES:
			EmitSamplesMode(bind_data, gs, output);
			return;
		case PcaMode::PCS:
			EmitPcsMode(bind_data, gs, output);
			return;
		case PcaMode::BOTH:
			EmitBothMode(bind_data, gs, output);
			return;
		}
	}

	// --- Not participating ---
	if (!ls.initialized || bind_data.effective_variants.empty() || ls.thread_id >= gs.thread_count) {
		output.SetCardinality(0);
		return;
	}

	// --- Wait for pass transition ---
	// If pass_generation hasn't advanced since our last participation,
	// the last thread is still doing merge/SVD. Return empty and wait.
	uint32_t gen = gs.pass_generation.load(std::memory_order_acquire);
	if (gen == ls.last_generation_seen) {
		output.SetCardinality(0);
		return;
	}
	ls.last_generation_seen = gen;

	// pass_generation starts at 0 and the first Scan entry sees
	// last_generation_seen == UINT32_MAX != 0, so it proceeds.
	// gen maps directly to the pass number.
	uint32_t pass = gen;
	if (pass >= gs.total_passes) {
		output.SetCardinality(0);
		return;
	}

	const uintptr_t *sample_include = nullptr;
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		sample_include = bind_data.sample_subset->SampleInclude();
	}

	// Register this thread as active in the current pass
	gs.pass_active_threads.fetch_add(1, std::memory_order_acq_rel);

	// Zero partial and scan variant blocks (parallel with other threads)
	std::fill(gs.thread_partials[ls.thread_id].begin(), gs.thread_partials[ls.thread_id].end(), 0.0);
	ScanVariantPass(bind_data, gs, ls, sample_include, pass);

	// --- Barrier ---
	uint32_t remaining = gs.pass_active_threads.fetch_sub(1, std::memory_order_acq_rel);
	if (remaining != 1) {
		// Not the last thread — return empty. DuckDB will re-call us and
		// we'll wait on pass_generation at the top until the last thread
		// finishes the merge and advances the generation.
		output.SetCardinality(0);
		return;
	}

	// --- Last thread: merge this pass, then loop through remaining passes ---
	// The last thread must not return 0 cardinality between passes because
	// DuckDB would interpret "all threads returned 0" as "scan complete."
	// Instead, the last thread processes the merge + next pass inline,
	// only returning when it either emits rows (algorithm done) or yields
	// to let other threads help with the next pass.
	while (true) {
		MergePass(bind_data, gs, pass);

		if (gs.algorithm_done.load(std::memory_order_relaxed)) {
			switch (bind_data.mode) {
			case PcaMode::SAMPLES:
				EmitSamplesMode(bind_data, gs, output);
				return;
			case PcaMode::PCS:
				EmitPcsMode(bind_data, gs, output);
				return;
			case PcaMode::BOTH:
				EmitBothMode(bind_data, gs, output);
				return;
			}
		}

		// Prepare next pass: reset block counter, advance generation to
		// release waiting threads. Order matters: reset state BEFORE
		// advancing generation so threads see clean state.
		pass++;
		gs.next_block_idx.store(0, std::memory_order_release);
		gs.pass_generation.fetch_add(1, std::memory_order_release);

		// Process the next pass inline. Register as active (matching
		// the normal Scan entry path) so the barrier works correctly
		// regardless of whether other threads have joined yet.
		gs.pass_active_threads.fetch_add(1, std::memory_order_acq_rel);
		std::fill(gs.thread_partials[ls.thread_id].begin(), gs.thread_partials[ls.thread_id].end(), 0.0);
		ScanVariantPass(bind_data, gs, ls, sample_include, pass);

		// Decrement and check if we're still last
		remaining = gs.pass_active_threads.fetch_sub(1, std::memory_order_acq_rel);

		if (remaining != 1) {
			// Other threads are still working on this pass. Return empty
			// and let DuckDB re-call us. On re-entry, we'll see
			// pass_generation == last_generation_seen and wait until
			// the actual last thread finishes and merges.
			ls.last_generation_seen = gs.pass_generation.load(std::memory_order_acquire);
			output.SetCardinality(0);
			return;
		}

		// We're still the last thread — loop back to merge
	}
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkPca(ExtensionLoader &loader) {
	TableFunction plink_pca("plink_pca", {LogicalType::VARCHAR}, PlinkPcaScan, PlinkPcaBind, PlinkPcaInitGlobal,
	                        PlinkPcaInitLocal);

	plink_pca.projection_pushdown = true;

	plink_pca.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_pca.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_pca.named_parameters["mode"] = LogicalType::VARCHAR;
	plink_pca.named_parameters["n_pcs"] = LogicalType::INTEGER;
	plink_pca.named_parameters["samples"] = LogicalType::ANY;
	plink_pca.named_parameters["region"] = LogicalType::VARCHAR;

	loader.RegisterFunction(plink_pca);
}

} // namespace duckdb
