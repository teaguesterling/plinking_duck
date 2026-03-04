#include "plink_score.hpp"
#include "plink_common.hpp"

#include "duckdb/parallel/task_scheduler.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <exception>
#include <mutex>
#include <thread>
#include <unordered_map>

namespace duckdb {

// ---------------------------------------------------------------------------
// Column indices
// ---------------------------------------------------------------------------

static constexpr idx_t COL_FID = 0;
static constexpr idx_t COL_IID = 1;
static constexpr idx_t COL_ALLELE_CT = 2;
static constexpr idx_t COL_DENOM = 3;
static constexpr idx_t COL_NAMED_ALLELE_DOSAGE_SUM = 4;
static constexpr idx_t COL_SCORE_SUM = 5;
static constexpr idx_t COL_SCORE_AVG = 6;

// ---------------------------------------------------------------------------
// ScoredVariant
// ---------------------------------------------------------------------------

struct ScoredVariant {
	uint32_t variant_idx; // index into .pgen file
	double weight;        // scoring weight
	bool flip;            // true if scored allele is REF (dosage = 2 - alt_dosage)
};

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkScoreBindData : public TableFunctionData {
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

	// Scored variants (built from weights parameter)
	vector<ScoredVariant> scored_variants;

	// Options
	bool center = false;
	bool no_mean_imputation = false;

	// Maps output index → original sample index (for sample metadata lookup)
	vector<uint32_t> sample_output_order;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkScoreGlobalState : public GlobalTableFunctionState {
	// Per-sample accumulators (filled during phase 1)
	vector<double> score_sums;
	vector<double> named_allele_dosage_sums;
	vector<uint32_t> allele_cts;

	// Phase 1 synchronization
	std::mutex scoring_mutex;
	std::atomic<bool> scoring_done {false};

	// Phase 2 emission
	std::atomic<uint32_t> next_sample_idx {0};
	uint32_t total_samples = 0;

	vector<column_t> column_ids;

	// DuckDB-configured thread count for internal parallelism
	uint32_t max_thread_count = 1;

	idx_t MaxThreads() const override {
		return 1;
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkScoreLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	// Buffers for PgrGetD
	AlignedBuffer genovec_buf;
	AlignedBuffer dosage_present_buf;
	AlignedBuffer dosage_main_buf;
	vector<double> dosage_doubles;

	bool initialized = false;

	~PlinkScoreLocalState() {
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

static unique_ptr<FunctionData> PlinkScoreBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkScoreBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters (first pass: file paths and options) ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "center") {
			bind_data->center = kv.second.GetValue<bool>();
		} else if (kv.first == "no_mean_imputation") {
			bind_data->no_mean_imputation = kv.second.GetValue<bool>();
		} else if (kv.first == "weights" || kv.first == "samples" || kv.first == "region") {
			// Handled below
		}
	}

	if (bind_data->center && bind_data->no_mean_imputation) {
		throw InvalidInputException("plink_score: center and no_mean_imputation cannot both be true");
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_score: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		if (bind_data->psam_path.empty()) {
			throw InvalidInputException("plink_score: cannot find .psam or .fam companion for '%s' "
			                            "(use psam := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
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
		throw IOException("plink_score: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
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
		throw IOException("plink_score: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadataIndex(context, bind_data->pvar_path, "plink_score");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_score: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (required for plink_score) ---
	bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);

	if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
		throw InvalidInputException("plink_score: sample count mismatch: .pgen has %u samples, "
		                            ".psam/.fam '%s' has %llu samples",
		                            bind_data->raw_sample_ct, bind_data->psam_path,
		                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
	}

	// --- Process samples parameter ---
	bind_data->effective_sample_ct = bind_data->raw_sample_ct;

	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		auto indices =
		    ResolveSampleIndices(samples_it->second, bind_data->raw_sample_ct, &bind_data->sample_info, "plink_score");

		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;

		// Build sample_output_order: sorted indices matching sample_include bit order
		auto sorted_indices = indices;
		std::sort(sorted_indices.begin(), sorted_indices.end());
		bind_data->sample_output_order = std::move(sorted_indices);
	} else {
		// No subsetting: identity mapping
		bind_data->sample_output_order.resize(bind_data->raw_sample_ct);
		for (uint32_t i = 0; i < bind_data->raw_sample_ct; i++) {
			bind_data->sample_output_order[i] = i;
		}
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range =
		    ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_score");
	}

	// --- Process weights parameter ---
	auto weights_it = input.named_parameters.find("weights");
	if (weights_it == input.named_parameters.end()) {
		throw InvalidInputException("plink_score: weights parameter is required");
	}

	auto &weights_val = weights_it->second;
	auto &weights_type = weights_val.type();

	// Determine variant range for positional mode
	uint32_t range_start = bind_data->variant_range.has_filter ? bind_data->variant_range.start_idx : 0;
	uint32_t range_end =
	    bind_data->variant_range.has_filter ? bind_data->variant_range.end_idx : bind_data->raw_variant_ct;
	uint32_t variant_count = range_end - range_start;

	if (weights_type.id() == LogicalTypeId::LIST) {
		auto &child_type = ListType::GetChildType(weights_type);
		auto &children = ListValue::GetChildren(weights_val);

		if (children.empty()) {
			throw InvalidInputException("plink_score: weights list is empty");
		}

		if (child_type.id() == LogicalTypeId::STRUCT) {
			// --- ID-keyed mode: LIST(STRUCT(id VARCHAR, allele VARCHAR, weight DOUBLE)) ---
			auto &struct_children = StructType::GetChildTypes(child_type);

			// Validate struct field names
			bool has_id = false, has_allele = false, has_weight = false;
			idx_t id_idx = 0, allele_idx = 0, weight_idx = 0;
			for (idx_t i = 0; i < struct_children.size(); i++) {
				if (struct_children[i].first == "id") {
					has_id = true;
					id_idx = i;
				} else if (struct_children[i].first == "allele") {
					has_allele = true;
					allele_idx = i;
				} else if (struct_children[i].first == "weight") {
					has_weight = true;
					weight_idx = i;
				}
			}

			if (!has_id || !has_allele || !has_weight) {
				throw InvalidInputException("plink_score: ID-keyed weights must be "
				                            "LIST(STRUCT(id VARCHAR, allele VARCHAR, weight DOUBLE))");
			}

			// Build variant ID → index map (respecting region filter)
			unordered_map<string, uint32_t> variant_id_map;
			for (uint32_t i = range_start; i < range_end; i++) {
				auto vid = bind_data->variants.GetId(i);
				if (!vid.empty()) {
					variant_id_map[vid] = i;
				}
			}

			uint32_t unmatched_id_count = 0;
			uint32_t unmatched_allele_count = 0;

			for (auto &entry : children) {
				auto &struct_vals = StructValue::GetChildren(entry);
				string id = struct_vals[id_idx].GetValue<string>();
				string allele = struct_vals[allele_idx].GetValue<string>();
				double weight = struct_vals[weight_idx].GetValue<double>();

				auto it = variant_id_map.find(id);
				if (it == variant_id_map.end()) {
					unmatched_id_count++;
					continue;
				}
				uint32_t vidx = it->second;

				// Determine allele orientation
				bool flip = false;
				if (allele == bind_data->variants.GetAlt(vidx)) {
					flip = false; // ALT allele: use dosage as-is
				} else if (allele == bind_data->variants.GetRef(vidx)) {
					flip = true; // REF allele: scored_dosage = 2 - alt_dosage
				} else {
					unmatched_allele_count++;
					continue;
				}

				if (weight != 0.0) {
					bind_data->scored_variants.push_back({vidx, weight, flip});
				}
			}

			// Sort by variant index for sequential .pgen access
			std::sort(bind_data->scored_variants.begin(), bind_data->scored_variants.end(),
			          [](const ScoredVariant &a, const ScoredVariant &b) { return a.variant_idx < b.variant_idx; });

		} else {
			// --- Positional mode: LIST(numeric) ---
			if (static_cast<uint32_t>(children.size()) != variant_count) {
				throw InvalidInputException("plink_score: weights list length (%llu) must match variant count (%u)",
				                            static_cast<unsigned long long>(children.size()), variant_count);
			}

			for (idx_t i = 0; i < children.size(); i++) {
				double w = children[i].GetValue<double>();
				if (w != 0.0) {
					bind_data->scored_variants.push_back({range_start + static_cast<uint32_t>(i), w, false});
				}
			}
		}
	} else {
		throw InvalidInputException("plink_score: weights must be a list (LIST(DOUBLE) for positional mode, "
		                            "or LIST(STRUCT(id, allele, weight)) for ID-keyed mode)");
	}

	// --- Register output columns ---
	names = {"FID", "IID", "ALLELE_CT", "DENOM", "NAMED_ALLELE_DOSAGE_SUM", "SCORE_SUM", "SCORE_AVG"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	                LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE};

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkScoreInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkScoreBindData>();
	auto state = make_uniq<PlinkScoreGlobalState>();

	state->total_samples = bind_data.effective_sample_ct;
	state->column_ids = input.column_ids;
	state->max_thread_count =
	    static_cast<uint32_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());

	// Initialize accumulators
	state->score_sums.resize(state->total_samples, 0.0);
	state->named_allele_dosage_sums.resize(state->total_samples, 0.0);
	state->allele_cts.resize(state->total_samples, 0);

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkScoreInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                               GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkScoreBindData>();
	auto state = make_uniq<PlinkScoreLocalState>();

	if (bind_data.scored_variants.empty()) {
		// No variants to score — skip pgenlib init
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
		throw IOException("plink_score: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_score: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_score: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	// Allocate buffers for PgrGetD
	uint32_t raw_sample_ct = bind_data.raw_sample_ct;
	uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(raw_sample_ct);
	state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));
	std::memset(state->genovec_buf.ptr, 0, genovec_word_ct * sizeof(uintptr_t));

	uintptr_t dosage_present_word_ct = plink2::BitCtToAlignedWordCt(raw_sample_ct);
	state->dosage_present_buf.Allocate(dosage_present_word_ct * sizeof(uintptr_t));
	std::memset(state->dosage_present_buf.ptr, 0, dosage_present_word_ct * sizeof(uintptr_t));

	state->dosage_main_buf.Allocate(raw_sample_ct * sizeof(uint16_t));
	std::memset(state->dosage_main_buf.ptr, 0, raw_sample_ct * sizeof(uint16_t));

	state->dosage_doubles.resize(bind_data.effective_sample_ct, 0.0);

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Scoring phase (called once, guarded by mutex)
// ---------------------------------------------------------------------------

static constexpr uint32_t SCORE_BATCH_SIZE = 16;

//! Per-thread scoring accumulators.
struct ScoreAccumulator {
	vector<double> score_sums;
	vector<double> dosage_sums;
	vector<uint32_t> allele_cts;

	void Init(uint32_t sample_ct) {
		score_sums.resize(sample_ct, 0.0);
		dosage_sums.resize(sample_ct, 0.0);
		allele_cts.resize(sample_ct, 0);
	}
};

//! Worker function for parallel scoring. Each worker claims batches from
//! scored_variants by atomic index, reads genotypes via its own PgenReader,
//! and accumulates into its thread-local ScoreAccumulator.
static void ScoreWorkerFunc(plink2::PgenReader *pgr, plink2::PgrSampleSubsetIndex pssi,
                            const uintptr_t *sample_include, uint32_t sample_ct, uintptr_t *genovec_buf,
                            uintptr_t *dosage_present_buf, uint16_t *dosage_main_buf, vector<double> &dosage_doubles,
                            const vector<ScoredVariant> &scored_variants, bool center, bool no_mean_imputation,
                            std::atomic<uint32_t> &work_counter, std::atomic<bool> &stop_requested,
                            ScoreAccumulator &accum, std::exception_ptr &exc_ptr) {
	uint32_t total_scored = static_cast<uint32_t>(scored_variants.size());
	try {
		while (!stop_requested.load(std::memory_order_relaxed)) {
			uint32_t batch_start = work_counter.fetch_add(SCORE_BATCH_SIZE);
			if (batch_start >= total_scored) {
				break;
			}
			uint32_t batch_end = std::min(batch_start + SCORE_BATCH_SIZE, total_scored);

			for (uint32_t si = batch_start; si < batch_end; si++) {
				auto &sv = scored_variants[si];
				uint32_t dosage_ct = 0;

				plink2::PglErr err = plink2::PgrGetD(sample_include, pssi, sample_ct, sv.variant_idx, pgr, genovec_buf,
				                                     dosage_present_buf, dosage_main_buf, &dosage_ct);
				if (err != plink2::kPglRetSuccess) {
					throw IOException("plink_score: PgrGetD failed for variant %u", sv.variant_idx);
				}

				plink2::Dosage16ToDoublesMinus9(genovec_buf, dosage_present_buf, dosage_main_buf, sample_ct, dosage_ct,
				                                dosage_doubles.data());

				// Per-variant statistics
				double sum_alt = 0.0;
				uint32_t non_missing_ct = 0;
				for (uint32_t s = 0; s < sample_ct; s++) {
					if (dosage_doubles[s] != -9.0) {
						sum_alt += dosage_doubles[s];
						non_missing_ct++;
					}
				}

				if (non_missing_ct == 0) {
					continue;
				}

				if (center) {
					double mean_alt = sum_alt / static_cast<double>(non_missing_ct);
					double freq = mean_alt / 2.0;
					double sd = std::sqrt(2.0 * freq * (1.0 - freq));
					if (sd == 0.0) {
						continue;
					}
					double mean_scored = sv.flip ? (2.0 - mean_alt) : mean_alt;

					for (uint32_t s = 0; s < sample_ct; s++) {
						if (dosage_doubles[s] == -9.0) {
							continue;
						}
						double scored_dosage = sv.flip ? (2.0 - dosage_doubles[s]) : dosage_doubles[s];
						double standardized = (scored_dosage - mean_scored) / sd;
						accum.score_sums[s] += sv.weight * standardized;
						accum.allele_cts[s] += 2;
					}
				} else if (no_mean_imputation) {
					for (uint32_t s = 0; s < sample_ct; s++) {
						if (dosage_doubles[s] == -9.0) {
							continue;
						}
						double scored_dosage = sv.flip ? (2.0 - dosage_doubles[s]) : dosage_doubles[s];
						accum.score_sums[s] += sv.weight * scored_dosage;
						accum.dosage_sums[s] += scored_dosage;
						accum.allele_cts[s] += 2;
					}
				} else {
					double mean_alt = sum_alt / static_cast<double>(non_missing_ct);
					for (uint32_t s = 0; s < sample_ct; s++) {
						double alt_dosage = (dosage_doubles[s] == -9.0) ? mean_alt : dosage_doubles[s];
						double scored_dosage = sv.flip ? (2.0 - alt_dosage) : alt_dosage;
						accum.score_sums[s] += sv.weight * scored_dosage;
						accum.dosage_sums[s] += scored_dosage;
						accum.allele_cts[s] += 2;
					}
				}
			}
		}
	} catch (...) {
		exc_ptr = std::current_exception();
		stop_requested.store(true, std::memory_order_relaxed);
	}
}

static void PerformScoring(const PlinkScoreBindData &bind_data, PlinkScoreGlobalState &gstate,
                           PlinkScoreLocalState &lstate) {
	uint32_t sample_ct = bind_data.effective_sample_ct;
	uint32_t scored_count = static_cast<uint32_t>(bind_data.scored_variants.size());

	const uintptr_t *sample_include = nullptr;
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		sample_include = bind_data.sample_subset->SampleInclude();
	}

	uint32_t db_threads = gstate.max_thread_count;
	uint32_t n_threads = (scored_count >= 100 && db_threads > 1)
	                         ? std::min({db_threads, scored_count / 16 + 1, static_cast<uint32_t>(16)})
	                         : 1;

	if (n_threads <= 1) {
		// Single-threaded path: accumulate directly into gstate
		ScoreAccumulator accum;
		accum.score_sums.assign(gstate.score_sums.begin(), gstate.score_sums.end());
		accum.dosage_sums.assign(gstate.named_allele_dosage_sums.begin(), gstate.named_allele_dosage_sums.end());
		accum.allele_cts.assign(gstate.allele_cts.begin(), gstate.allele_cts.end());

		std::atomic<uint32_t> work_counter(0);
		std::atomic<bool> stop_requested(false);
		std::exception_ptr exc_ptr = nullptr;

		ScoreWorkerFunc(&lstate.pgr, lstate.pssi, sample_include, sample_ct, lstate.genovec_buf.As<uintptr_t>(),
		                lstate.dosage_present_buf.As<uintptr_t>(), lstate.dosage_main_buf.As<uint16_t>(),
		                lstate.dosage_doubles, bind_data.scored_variants, bind_data.center,
		                bind_data.no_mean_imputation, work_counter, stop_requested, accum, exc_ptr);

		if (exc_ptr) {
			std::rethrow_exception(exc_ptr);
		}

		gstate.score_sums = std::move(accum.score_sums);
		gstate.named_allele_dosage_sums = std::move(accum.dosage_sums);
		gstate.allele_cts = std::move(accum.allele_cts);
	} else {
		// Parallel path
		std::atomic<uint32_t> work_counter(0);
		std::atomic<bool> stop_requested(false);

		vector<ScoreAccumulator> accums(n_threads);
		for (uint32_t t = 0; t < n_threads; t++) {
			accums[t].Init(sample_ct);
		}

		// Per-thread pgenlib contexts and buffers
		uint32_t raw_sample_ct = bind_data.raw_sample_ct;
		uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(raw_sample_ct);
		uintptr_t dosage_present_word_ct = plink2::BitCtToAlignedWordCt(raw_sample_ct);

		vector<PgenThreadContext> thread_ctxs(n_threads - 1);
		vector<AlignedBuffer> thread_genovec_bufs(n_threads - 1);
		vector<AlignedBuffer> thread_dosage_present_bufs(n_threads - 1);
		vector<AlignedBuffer> thread_dosage_main_bufs(n_threads - 1);
		vector<vector<double>> thread_dosage_doubles(n_threads - 1);
		vector<std::exception_ptr> exc_ptrs(n_threads, nullptr);

		for (uint32_t t = 0; t < n_threads - 1; t++) {
			InitPgenThreadContext(thread_ctxs[t], bind_data.pgen_path, bind_data.raw_variant_ct, raw_sample_ct,
			                      bind_data.has_sample_subset ? bind_data.sample_subset.get() : nullptr);

			thread_genovec_bufs[t].Allocate(genovec_word_ct * sizeof(uintptr_t));
			std::memset(thread_genovec_bufs[t].ptr, 0, genovec_word_ct * sizeof(uintptr_t));

			thread_dosage_present_bufs[t].Allocate(dosage_present_word_ct * sizeof(uintptr_t));
			std::memset(thread_dosage_present_bufs[t].ptr, 0, dosage_present_word_ct * sizeof(uintptr_t));

			thread_dosage_main_bufs[t].Allocate(raw_sample_ct * sizeof(uint16_t));
			std::memset(thread_dosage_main_bufs[t].ptr, 0, raw_sample_ct * sizeof(uint16_t));

			thread_dosage_doubles[t].resize(sample_ct, 0.0);
		}

		// Launch N-1 background threads
		vector<std::thread> threads;
		threads.reserve(n_threads - 1);
		for (uint32_t t = 0; t < n_threads - 1; t++) {
			threads.emplace_back(ScoreWorkerFunc, &thread_ctxs[t].pgr, thread_ctxs[t].pssi, sample_include, sample_ct,
			                     thread_genovec_bufs[t].As<uintptr_t>(),
			                     thread_dosage_present_bufs[t].As<uintptr_t>(),
			                     thread_dosage_main_bufs[t].As<uint16_t>(), std::ref(thread_dosage_doubles[t]),
			                     std::cref(bind_data.scored_variants), bind_data.center, bind_data.no_mean_imputation,
			                     std::ref(work_counter), std::ref(stop_requested), std::ref(accums[t + 1]),
			                     std::ref(exc_ptrs[t + 1]));
		}

		// Main thread as worker 0
		ScoreWorkerFunc(&lstate.pgr, lstate.pssi, sample_include, sample_ct, lstate.genovec_buf.As<uintptr_t>(),
		                lstate.dosage_present_buf.As<uintptr_t>(), lstate.dosage_main_buf.As<uint16_t>(),
		                lstate.dosage_doubles, bind_data.scored_variants, bind_data.center,
		                bind_data.no_mean_imputation, work_counter, stop_requested, accums[0], exc_ptrs[0]);

		// Join all threads
		for (auto &t : threads) {
			t.join();
		}

		// Check for exceptions
		for (auto &eptr : exc_ptrs) {
			if (eptr) {
				std::rethrow_exception(eptr);
			}
		}

		// Merge accumulators into gstate
		for (uint32_t t = 0; t < n_threads; t++) {
			for (uint32_t s = 0; s < sample_ct; s++) {
				gstate.score_sums[s] += accums[t].score_sums[s];
				gstate.named_allele_dosage_sums[s] += accums[t].dosage_sums[s];
				gstate.allele_cts[s] += accums[t].allele_cts[s];
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void PlinkScoreScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkScoreBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkScoreGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkScoreLocalState>();

	// Phase 1: Score all variants (once, guarded by mutex)
	if (!gstate.scoring_done.load(std::memory_order_acquire)) {
		std::lock_guard<std::mutex> lock(gstate.scoring_mutex);
		if (!gstate.scoring_done.load(std::memory_order_relaxed)) {
			if (lstate.initialized && !bind_data.scored_variants.empty()) {
				PerformScoring(bind_data, gstate, lstate);
			}
			gstate.scoring_done.store(true, std::memory_order_release);
		}
	}

	// Phase 2: Emit one row per sample
	auto &column_ids = gstate.column_ids;
	uint32_t total_samples = gstate.total_samples;
	bool has_fid = !bind_data.sample_info.fids.empty();

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t sidx = gstate.next_sample_idx.fetch_add(1);
		if (sidx >= total_samples) {
			break;
		}

		// Map output index to original sample index
		uint32_t orig_idx = bind_data.sample_output_order[sidx];

		uint32_t allele_ct = gstate.allele_cts[sidx];
		double score_sum = gstate.score_sums[sidx];
		double dosage_sum = gstate.named_allele_dosage_sums[sidx];
		double score_avg = (allele_ct > 0) ? score_sum / static_cast<double>(allele_ct) : 0.0;

		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];
			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];

			switch (file_col) {
			case COL_FID: {
				if (has_fid) {
					FlatVector::GetData<string_t>(vec)[rows_emitted] =
					    StringVector::AddString(vec, bind_data.sample_info.fids[orig_idx]);
				} else {
					FlatVector::SetNull(vec, rows_emitted, true);
				}
				break;
			}
			case COL_IID: {
				FlatVector::GetData<string_t>(vec)[rows_emitted] =
				    StringVector::AddString(vec, bind_data.sample_info.iids[orig_idx]);
				break;
			}
			case COL_ALLELE_CT: {
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(allele_ct);
				break;
			}
			case COL_DENOM: {
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(allele_ct);
				break;
			}
			case COL_NAMED_ALLELE_DOSAGE_SUM: {
				FlatVector::GetData<double>(vec)[rows_emitted] = dosage_sum;
				break;
			}
			case COL_SCORE_SUM: {
				FlatVector::GetData<double>(vec)[rows_emitted] = score_sum;
				break;
			}
			case COL_SCORE_AVG: {
				FlatVector::GetData<double>(vec)[rows_emitted] = score_avg;
				break;
			}
			default:
				break;
			}
		}

		rows_emitted++;
	}

	output.SetCardinality(rows_emitted);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkScore(ExtensionLoader &loader) {
	TableFunction plink_score("plink_score", {LogicalType::VARCHAR}, PlinkScoreScan, PlinkScoreBind,
	                          PlinkScoreInitGlobal, PlinkScoreInitLocal);

	plink_score.projection_pushdown = true;

	plink_score.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_score.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_score.named_parameters["weights"] = LogicalType::ANY;
	plink_score.named_parameters["samples"] = LogicalType::ANY;
	plink_score.named_parameters["region"] = LogicalType::VARCHAR;
	plink_score.named_parameters["center"] = LogicalType::BOOLEAN;
	plink_score.named_parameters["no_mean_imputation"] = LogicalType::BOOLEAN;

	loader.RegisterFunction(plink_score);
}

} // namespace duckdb
