#include "plink_ld.hpp"
#include "plink_common.hpp"

#include <atomic>
#include <cmath>

namespace duckdb {

// ---------------------------------------------------------------------------
// Column indices
// ---------------------------------------------------------------------------

static constexpr idx_t COL_CHROM_A = 0;
static constexpr idx_t COL_POS_A = 1;
static constexpr idx_t COL_ID_A = 2;
static constexpr idx_t COL_CHROM_B = 3;
static constexpr idx_t COL_POS_B = 4;
static constexpr idx_t COL_ID_B = 5;
static constexpr idx_t COL_R2 = 6;
static constexpr idx_t COL_D_PRIME = 7;
static constexpr idx_t COL_OBS_CT = 8;

// ---------------------------------------------------------------------------
// LD computation
// ---------------------------------------------------------------------------

enum class LdMode : uint8_t { PAIRWISE, WINDOWED };

struct LdResult {
	double r2;
	double d_prime;
	uint32_t obs_ct;
	bool is_valid; // false if monomorphic, < 2 obs, etc.
};

//! Compute LD statistics from two packed 2-bit genotype arrays.
//! Genotype encoding: 0=hom_ref, 1=het, 2=hom_alt, 3=missing.
static LdResult ComputeLdStats(const uintptr_t *genovec_a, const uintptr_t *genovec_b, uint32_t sample_ct) {
	double sum_a = 0, sum_b = 0, sum_ab = 0, sum_a2 = 0, sum_b2 = 0;
	uint32_t n = 0;

	uint32_t word_ct = plink2::DivUp(sample_ct, plink2::kBitsPerWordD2);
	for (uint32_t widx = 0; widx < word_ct; widx++) {
		uintptr_t word_a = genovec_a[widx];
		uintptr_t word_b = genovec_b[widx];

		uint32_t samples_remaining = sample_ct - widx * plink2::kBitsPerWordD2;
		uint32_t samples_in_word = std::min(samples_remaining, static_cast<uint32_t>(plink2::kBitsPerWordD2));

		for (uint32_t sidx = 0; sidx < samples_in_word; sidx++) {
			uint32_t geno_a = word_a & 3;
			uint32_t geno_b = word_b & 3;
			word_a >>= 2;
			word_b >>= 2;

			if (geno_a == 3 || geno_b == 3) {
				continue;
			}

			double ga = static_cast<double>(geno_a);
			double gb = static_cast<double>(geno_b);
			sum_a += ga;
			sum_b += gb;
			sum_ab += ga * gb;
			sum_a2 += ga * ga;
			sum_b2 += gb * gb;
			n++;
		}
	}

	LdResult result;
	result.obs_ct = n;
	result.is_valid = false;
	result.r2 = 0;
	result.d_prime = 0;

	if (n < 2) {
		return result;
	}

	double dn = static_cast<double>(n);
	double mean_a = sum_a / dn;
	double mean_b = sum_b / dn;
	double cov_ab = sum_ab / dn - mean_a * mean_b;
	double var_a = sum_a2 / dn - mean_a * mean_a;
	double var_b = sum_b2 / dn - mean_b * mean_b;

	// Monomorphic variant — correlation undefined
	if (var_a < 1e-15 || var_b < 1e-15) {
		return result;
	}

	result.is_valid = true;
	result.r2 = (cov_ab * cov_ab) / (var_a * var_b);

	// D' via composite LD estimator (Weir 1979):
	//   D = cov(gA, gB) / 4
	//   D' = D / D_max where D_max depends on sign of D
	// Note: this estimator uses genotype-level (not haplotype-level) statistics,
	// so D' can exceed 1.0 when samples deviate from Hardy-Weinberg equilibrium.
	double D = cov_ab / 4.0;
	double p_a = sum_a / (2.0 * dn);
	double p_b = sum_b / (2.0 * dn);

	double D_max;
	if (D >= 0) {
		D_max = std::min(p_a * (1.0 - p_b), (1.0 - p_a) * p_b);
	} else {
		D_max = std::max(-p_a * p_b, -(1.0 - p_a) * (1.0 - p_b));
	}

	if (std::abs(D_max) < 1e-15) {
		result.d_prime = 0.0;
	} else {
		// D/D_max is always non-negative with this formula
		result.d_prime = D / D_max;
	}

	return result;
}

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkLdBindData : public TableFunctionData {
	string pgen_path;
	string pvar_path;
	string psam_path;

	VariantMetadataIndex variants;
	SampleInfo sample_info;
	bool has_sample_info = false;

	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Sample subsetting
	bool has_sample_subset = false;
	unique_ptr<SampleSubset> sample_subset;
	uint32_t effective_sample_ct = 0;

	// Region filtering
	VariantRange variant_range;

	// Mode
	LdMode mode = LdMode::WINDOWED;

	// Pairwise mode
	uint32_t pairwise_vidx_a = 0;
	uint32_t pairwise_vidx_b = 0;

	// Windowed mode parameters
	int64_t window_bp = 1000000; // window_kb * 1000
	double r2_threshold = 0.2;
	bool inter_chr = false;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkLdGlobalState : public GlobalTableFunctionState {
	LdMode mode;
	uint32_t start_variant_idx = 0;
	uint32_t end_variant_idx = 0;

	// Pairwise mode
	std::atomic<bool> pair_emitted {false};

	// Windowed mode
	std::atomic<uint32_t> next_anchor_idx {0};

	idx_t MaxThreads() const override {
		if (mode == LdMode::PAIRWISE) {
			return 1;
		}
		uint32_t range = end_variant_idx - start_variant_idx;
		return std::min<idx_t>(range / 50 + 1, 16);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkLdLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	AlignedBuffer genovec_a_buf; // anchor variant genotypes
	AlignedBuffer genovec_b_buf; // partner variant genotypes

	// Windowed mode: state preservation across scan calls
	bool in_window = false;
	uint32_t anchor_idx = 0;
	uint32_t next_j = 0;

	bool initialized = false;

	~PlinkLdLocalState() {
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

static unique_ptr<FunctionData> PlinkLdBind(ClientContext &context, TableFunctionBindInput &input,
                                            vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkLdBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Collect named parameters ---
	string variant1_id, variant2_id;

	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "variant1") {
			variant1_id = kv.second.GetValue<string>();
		} else if (kv.first == "variant2") {
			variant2_id = kv.second.GetValue<string>();
		} else if (kv.first == "window_kb") {
			auto kb = kv.second.GetValue<int64_t>();
			if (kb < 0) {
				throw InvalidInputException("plink_ld: window_kb must be non-negative");
			}
			bind_data->window_bp = kb * 1000;
		} else if (kv.first == "r2_threshold") {
			bind_data->r2_threshold = kv.second.GetValue<double>();
			if (bind_data->r2_threshold < 0.0 || bind_data->r2_threshold > 1.0) {
				throw InvalidInputException("plink_ld: r2_threshold must be between 0.0 and 1.0");
			}
		} else if (kv.first == "inter_chr") {
			bind_data->inter_chr = kv.second.GetValue<bool>();
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled after pgenlib init
		}
	}

	// --- Determine mode ---
	if (!variant1_id.empty() && !variant2_id.empty()) {
		bind_data->mode = LdMode::PAIRWISE;
	} else if (!variant1_id.empty() || !variant2_id.empty()) {
		throw InvalidInputException("plink_ld: both variant1 and variant2 must be specified for pairwise mode");
	} else {
		bind_data->mode = LdMode::WINDOWED;
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_ld: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		// .psam is optional for plink_ld — only needed if samples parameter uses VARCHAR IDs
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
		throw IOException("plink_ld: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
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
		throw IOException("plink_ld: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadataIndex(context, bind_data->pvar_path, "plink_ld");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_ld: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (optional) ---
	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);
		bind_data->has_sample_info = true;

		if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_ld: sample count mismatch: .pgen has %u samples, "
			                            ".psam/.fam '%s' has %llu samples",
			                            bind_data->raw_sample_ct, bind_data->psam_path,
			                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
		}
	}

	// --- Process samples parameter ---
	bind_data->effective_sample_ct = bind_data->raw_sample_ct;

	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		auto indices = ResolveSampleIndices(samples_it->second, bind_data->raw_sample_ct,
		                                    bind_data->has_sample_info ? &bind_data->sample_info : nullptr, "plink_ld");

		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range = ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_ld");
	}

	// --- Resolve pairwise variant indices ---
	if (bind_data->mode == LdMode::PAIRWISE) {
		// Look up variant1 by ID
		bool found_a = false, found_b = false;
		for (uint32_t i = 0; i < static_cast<uint32_t>(bind_data->variants.variant_ct); i++) {
			auto vid = bind_data->variants.GetId(i);
			if (vid == variant1_id) {
				bind_data->pairwise_vidx_a = i;
				found_a = true;
			}
			if (vid == variant2_id) {
				bind_data->pairwise_vidx_b = i;
				found_b = true;
			}
			if (found_a && found_b) {
				break;
			}
		}
		if (!found_a) {
			throw InvalidInputException("plink_ld: variant '%s' not found in .pvar", variant1_id);
		}
		if (!found_b) {
			throw InvalidInputException("plink_ld: variant '%s' not found in .pvar", variant2_id);
		}
	}

	// --- Register output columns ---
	names = {"CHROM_A", "POS_A", "ID_A", "CHROM_B", "POS_B", "ID_B", "R2", "D_PRIME", "OBS_CT"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
	                LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::INTEGER};

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkLdInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkLdBindData>();
	auto state = make_uniq<PlinkLdGlobalState>();

	state->mode = bind_data.mode;

	if (bind_data.variant_range.has_filter) {
		state->start_variant_idx = bind_data.variant_range.start_idx;
		state->end_variant_idx = bind_data.variant_range.end_idx;
	} else {
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.raw_variant_ct;
	}

	if (state->mode == LdMode::WINDOWED) {
		state->next_anchor_idx.store(state->start_variant_idx);
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkLdInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                            GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkLdBindData>();
	auto state = make_uniq<PlinkLdLocalState>();

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
		throw IOException("plink_ld: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_ld: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_ld: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	// Allocate two genovec buffers (anchor + partner)
	uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(bind_data.effective_sample_ct);
	uintptr_t genovec_bytes = genovec_word_ct * sizeof(uintptr_t);
	state->genovec_a_buf.Allocate(genovec_bytes);
	state->genovec_b_buf.Allocate(genovec_bytes);
	std::memset(state->genovec_a_buf.ptr, 0, genovec_bytes);
	std::memset(state->genovec_b_buf.ptr, 0, genovec_bytes);

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Helpers: emit an output row and read genotypes
// ---------------------------------------------------------------------------

static void EmitRow(DataChunk &output, idx_t row_idx, const PlinkLdBindData &bind_data, uint32_t vidx_a,
                    uint32_t vidx_b, const LdResult &result) {
	auto &variants = bind_data.variants;

	// CHROM_A
	FlatVector::GetData<string_t>(output.data[COL_CHROM_A])[row_idx] =
	    StringVector::AddString(output.data[COL_CHROM_A], variants.GetChrom(vidx_a));

	// POS_A
	FlatVector::GetData<int32_t>(output.data[COL_POS_A])[row_idx] = variants.GetPos(vidx_a);

	// ID_A
	auto id_a = variants.GetId(vidx_a);
	if (id_a.empty()) {
		FlatVector::SetNull(output.data[COL_ID_A], row_idx, true);
	} else {
		FlatVector::GetData<string_t>(output.data[COL_ID_A])[row_idx] =
		    StringVector::AddString(output.data[COL_ID_A], id_a);
	}

	// CHROM_B
	FlatVector::GetData<string_t>(output.data[COL_CHROM_B])[row_idx] =
	    StringVector::AddString(output.data[COL_CHROM_B], variants.GetChrom(vidx_b));

	// POS_B
	FlatVector::GetData<int32_t>(output.data[COL_POS_B])[row_idx] = variants.GetPos(vidx_b);

	// ID_B
	auto id_b = variants.GetId(vidx_b);
	if (id_b.empty()) {
		FlatVector::SetNull(output.data[COL_ID_B], row_idx, true);
	} else {
		FlatVector::GetData<string_t>(output.data[COL_ID_B])[row_idx] =
		    StringVector::AddString(output.data[COL_ID_B], id_b);
	}

	// R2, D_PRIME, OBS_CT
	if (result.is_valid) {
		FlatVector::GetData<double>(output.data[COL_R2])[row_idx] = result.r2;
		FlatVector::GetData<double>(output.data[COL_D_PRIME])[row_idx] = result.d_prime;
	} else {
		FlatVector::SetNull(output.data[COL_R2], row_idx, true);
		FlatVector::SetNull(output.data[COL_D_PRIME], row_idx, true);
	}
	FlatVector::GetData<int32_t>(output.data[COL_OBS_CT])[row_idx] = static_cast<int32_t>(result.obs_ct);
}

static void ReadGenovec(PlinkLdLocalState &lstate, const PlinkLdBindData &bind_data, uint32_t vidx,
                        uintptr_t *genovec_out) {
	const uintptr_t *sample_include = nullptr;
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		sample_include = bind_data.sample_subset->SampleInclude();
	}

	plink2::PglErr err =
	    plink2::PgrGet(sample_include, lstate.pssi, bind_data.effective_sample_ct, vidx, &lstate.pgr, genovec_out);

	if (err != plink2::kPglRetSuccess) {
		throw IOException("plink_ld: PgrGet failed for variant %u", vidx);
	}
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void PlinkLdScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkLdBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkLdGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkLdLocalState>();

	if (!lstate.initialized) {
		output.SetCardinality(0);
		return;
	}

	uint32_t sample_ct = bind_data.effective_sample_ct;
	auto *genovec_a = lstate.genovec_a_buf.As<uintptr_t>();
	auto *genovec_b = lstate.genovec_b_buf.As<uintptr_t>();

	idx_t rows_emitted = 0;

	if (gstate.mode == LdMode::PAIRWISE) {
		// --- Pairwise mode: emit a single pair ---
		if (gstate.pair_emitted.exchange(true)) {
			output.SetCardinality(0);
			return;
		}

		uint32_t vidx_a = bind_data.pairwise_vidx_a;
		uint32_t vidx_b = bind_data.pairwise_vidx_b;

		ReadGenovec(lstate, bind_data, vidx_a, genovec_a);
		if (vidx_a == vidx_b) {
			// Self-LD: use same buffer for both
			auto result = ComputeLdStats(genovec_a, genovec_a, sample_ct);
			EmitRow(output, 0, bind_data, vidx_a, vidx_b, result);
		} else {
			ReadGenovec(lstate, bind_data, vidx_b, genovec_b);
			auto result = ComputeLdStats(genovec_a, genovec_b, sample_ct);
			EmitRow(output, 0, bind_data, vidx_a, vidx_b, result);
		}
		output.SetCardinality(1);
		return;
	}

	// --- Windowed mode ---
	uint32_t end_idx = gstate.end_variant_idx;
	auto &variants = bind_data.variants;

	// Cache anchor chrom/pos for the inner loop to avoid repeated parsing
	string anchor_chrom;
	int32_t anchor_pos = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		if (lstate.in_window) {
			// Resume scanning partners for current anchor
			uint32_t ai = lstate.anchor_idx;
			uint32_t j = lstate.next_j;

			if (anchor_chrom.empty()) {
				anchor_chrom = variants.GetChrom(ai);
				anchor_pos = variants.GetPos(ai);
			}

			while (j < end_idx) {
				auto j_chrom = variants.GetChrom(j);
				bool same_chrom = (j_chrom == anchor_chrom);

				if (same_chrom) {
					int64_t dist =
					    static_cast<int64_t>(variants.GetPos(j)) - static_cast<int64_t>(anchor_pos);
					if (dist > bind_data.window_bp) {
						if (!bind_data.inter_chr) {
							break; // Past window, no cross-chrom needed
						}
						// Skip remaining same-chrom variants beyond window
						while (j < end_idx && variants.GetChrom(j) == anchor_chrom) {
							j++;
						}
						continue; // Check next chromosome
					}
				} else {
					if (!bind_data.inter_chr) {
						break; // Past chromosome boundary
					}
					// Inter-chr: no distance filter
				}

				ReadGenovec(lstate, bind_data, j, genovec_b);
				auto result = ComputeLdStats(genovec_a, genovec_b, sample_ct);

				if (result.is_valid && result.r2 >= bind_data.r2_threshold) {
					EmitRow(output, rows_emitted, bind_data, ai, j, result);
					rows_emitted++;

					if (rows_emitted >= STANDARD_VECTOR_SIZE) {
						lstate.next_j = j + 1;
						// genovec_a still has anchor data
						goto done;
					}
				}

				j++;
			}

			// Anchor's window exhausted
			lstate.in_window = false;
		}

		// Claim a new anchor variant
		uint32_t anchor_idx = gstate.next_anchor_idx.fetch_add(1);
		if (anchor_idx >= end_idx) {
			break;
		}

		// Load anchor genotypes
		ReadGenovec(lstate, bind_data, anchor_idx, genovec_a);

		lstate.anchor_idx = anchor_idx;
		lstate.next_j = anchor_idx + 1;
		lstate.in_window = true;
		anchor_chrom = variants.GetChrom(anchor_idx);
		anchor_pos = variants.GetPos(anchor_idx);
	}

done:
	output.SetCardinality(rows_emitted);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkLd(ExtensionLoader &loader) {
	TableFunction plink_ld("plink_ld", {LogicalType::VARCHAR}, PlinkLdScan, PlinkLdBind, PlinkLdInitGlobal,
	                       PlinkLdInitLocal);

	plink_ld.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_ld.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_ld.named_parameters["variant1"] = LogicalType::VARCHAR;
	plink_ld.named_parameters["variant2"] = LogicalType::VARCHAR;
	plink_ld.named_parameters["window_kb"] = LogicalType::INTEGER;
	plink_ld.named_parameters["r2_threshold"] = LogicalType::DOUBLE;
	plink_ld.named_parameters["region"] = LogicalType::VARCHAR;
	plink_ld.named_parameters["samples"] = LogicalType::ANY;
	plink_ld.named_parameters["inter_chr"] = LogicalType::BOOLEAN;

	loader.RegisterFunction(plink_ld);
}

} // namespace duckdb
