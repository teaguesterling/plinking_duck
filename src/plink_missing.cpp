#include "plink_missing.hpp"
#include "plink_common.hpp"

#include <atomic>

namespace duckdb {

// ---------------------------------------------------------------------------
// Column indices — variant mode
// ---------------------------------------------------------------------------

// CHROM(0) POS(1) ID(2) REF(3) ALT(4) MISSING_CT(5) OBS_CT(6) F_MISS(7)
static constexpr idx_t VCOL_CHROM = 0;
static constexpr idx_t VCOL_POS = 1;
static constexpr idx_t VCOL_ID = 2;
static constexpr idx_t VCOL_REF = 3;
static constexpr idx_t VCOL_ALT = 4;
static constexpr idx_t VCOL_MISSING_CT = 5;
static constexpr idx_t VCOL_OBS_CT = 6;
static constexpr idx_t VCOL_F_MISS = 7;

// ---------------------------------------------------------------------------
// Column indices — sample mode
// ---------------------------------------------------------------------------

// FID(0) IID(1) MISSING_CT(2) OBS_CT(3) F_MISS(4)
static constexpr idx_t SCOL_FID = 0;
static constexpr idx_t SCOL_IID = 1;
static constexpr idx_t SCOL_MISSING_CT = 2;
static constexpr idx_t SCOL_OBS_CT = 3;
static constexpr idx_t SCOL_F_MISS = 4;

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkMissingBindData : public TableFunctionData {
	string pgen_path;
	string pvar_path;
	string psam_path;

	VariantMetadata variants;
	SampleInfo sample_info;
	bool has_sample_info = false;

	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Sample subsetting
	bool has_sample_subset = false;
	unique_ptr<SampleSubset> sample_subset;
	uint32_t effective_sample_ct = 0;

	// For sample mode with subsetting: maps subset position → original sample index
	vector<uint32_t> subset_original_indices;

	// Region filtering
	VariantRange variant_range;

	// Mode
	bool sample_mode = false;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkMissingGlobalState : public GlobalTableFunctionState {
	// Variant scanning range (used by both modes)
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t start_variant_idx = 0;
	uint32_t end_variant_idx = 0;

	// Projection pushdown
	vector<column_t> column_ids;
	bool need_missingness = false;

	// Mode
	bool sample_mode = false;

	// Sample mode: per-sample accumulation (currently single-threaded;
	// atomic for safety if MaxThreads is increased later)
	vector<uint32_t> sample_missing_counts;
	std::atomic<bool> variant_scan_done {false};
	std::atomic<uint32_t> next_sample_idx {0};
	uint32_t total_variant_ct = 0;

	idx_t MaxThreads() const override {
		if (sample_mode) {
			return 1;
		}
		uint32_t range = end_variant_idx - start_variant_idx;
		return std::min<idx_t>(range / 500 + 1, 16);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkMissingLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	// Buffers for PgrGetMissingness
	AlignedBuffer missingness_buf;
	AlignedBuffer genovec_buf;

	bool initialized = false;

	~PlinkMissingLocalState() {
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

static unique_ptr<FunctionData> PlinkMissingBind(ClientContext &context, TableFunctionBindInput &input,
                                                 vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkMissingBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "mode") {
			auto mode_str = kv.second.GetValue<string>();
			if (mode_str == "variant") {
				bind_data->sample_mode = false;
			} else if (mode_str == "sample") {
				bind_data->sample_mode = true;
			} else {
				throw InvalidInputException("plink_missing: mode must be 'variant' or 'sample', got '%s'", mode_str);
			}
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled after pgenlib init
		}
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_missing: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
	}

	// Sample mode requires .psam for IID output
	if (bind_data->sample_mode && bind_data->psam_path.empty()) {
		throw InvalidInputException("plink_missing: sample mode requires a .psam or .fam file "
		                            "(use psam := 'path' to specify explicitly)");
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
		throw IOException("plink_missing: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;

	// Phase 2 (validate file; per-thread readers re-init later)
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
		throw IOException("plink_missing: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "plink_missing");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_missing: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (optional for variant mode, required for sample mode) ---
	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);
		bind_data->has_sample_info = true;

		if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_missing: sample count mismatch: .pgen has %u samples, "
			                            ".psam/.fam '%s' has %llu samples",
			                            bind_data->raw_sample_ct, bind_data->psam_path,
			                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
		}
	}

	// --- Process samples parameter ---
	bind_data->effective_sample_ct = bind_data->raw_sample_ct;

	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		auto indices =
		    ResolveSampleIndices(samples_it->second, bind_data->raw_sample_ct,
		                         bind_data->has_sample_info ? &bind_data->sample_info : nullptr, "plink_missing");

		// Store sorted indices for sample mode output mapping
		bind_data->subset_original_indices = indices;
		std::sort(bind_data->subset_original_indices.begin(), bind_data->subset_original_indices.end());

		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range =
		    ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_missing");
	}

	// --- Register output columns based on mode ---
	if (bind_data->sample_mode) {
		names = {"FID", "IID", "MISSING_CT", "OBS_CT", "F_MISS"};
		return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
		                LogicalType::DOUBLE};
	} else {
		names = {"CHROM", "POS", "ID", "REF", "ALT", "MISSING_CT", "OBS_CT", "F_MISS"};
		return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
		                LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::DOUBLE};
	}

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkMissingInitGlobal(ClientContext &context,
                                                                   TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkMissingBindData>();
	auto state = make_uniq<PlinkMissingGlobalState>();

	state->sample_mode = bind_data.sample_mode;

	if (bind_data.variant_range.has_filter) {
		state->start_variant_idx = bind_data.variant_range.start_idx;
		state->end_variant_idx = bind_data.variant_range.end_idx;
	} else {
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.raw_variant_ct;
	}

	state->next_variant_idx.store(state->start_variant_idx);
	state->column_ids = input.column_ids;

	// Check projection pushdown: skip PgrGetMissingness if no missingness columns projected
	state->need_missingness = false;
	if (bind_data.sample_mode) {
		for (auto col_id : input.column_ids) {
			if (col_id != COLUMN_IDENTIFIER_ROW_ID && col_id >= SCOL_MISSING_CT) {
				state->need_missingness = true;
				break;
			}
		}
	} else {
		for (auto col_id : input.column_ids) {
			if (col_id != COLUMN_IDENTIFIER_ROW_ID && col_id >= VCOL_MISSING_CT) {
				state->need_missingness = true;
				break;
			}
		}
	}

	// Sample mode: pre-allocate accumulation array
	if (bind_data.sample_mode) {
		state->total_variant_ct = state->end_variant_idx - state->start_variant_idx;
		if (state->need_missingness) {
			state->sample_missing_counts.resize(bind_data.effective_sample_ct, 0);
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader + missingness buffers)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkMissingInitLocal(ExecutionContext &context,
                                                                 TableFunctionInitInput &input,
                                                                 GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkMissingBindData>();
	auto &gstate = global_state->Cast<PlinkMissingGlobalState>();
	auto state = make_uniq<PlinkMissingLocalState>();

	if (!gstate.need_missingness) {
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
		throw IOException("plink_missing: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_missing: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_missing: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	// Allocate missingness bitarray buffer (one bit per effective sample)
	uintptr_t missingness_alloc_ct = plink2::BitCtToAlignedWordCt(bind_data.effective_sample_ct);
	state->missingness_buf.Allocate(missingness_alloc_ct * sizeof(uintptr_t));

	// Allocate genovec scratch buffer (2 bits per raw sample, for pgenlib internal use)
	uintptr_t genovec_alloc_ct = plink2::NypCtToAlignedWordCt(bind_data.raw_sample_ct);
	state->genovec_buf.Allocate(genovec_alloc_ct * sizeof(uintptr_t));

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Scan — variant mode (parallel, same pattern as plink_freq)
// ---------------------------------------------------------------------------

static constexpr uint32_t MISSING_BATCH_SIZE = 128;

static void PlinkMissingScanVariant(const PlinkMissingBindData &bind_data, PlinkMissingGlobalState &gstate,
                                    PlinkMissingLocalState &lstate, DataChunk &output) {
	auto &column_ids = gstate.column_ids;
	uint32_t end_idx = gstate.end_variant_idx;
	uint32_t sample_ct = bind_data.effective_sample_ct;

	const uintptr_t *sample_include = nullptr;
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		sample_include = bind_data.sample_subset->SampleInclude();
	}

	uintptr_t word_ct = plink2::BitCtToWordCt(sample_ct);

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(MISSING_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
		if (batch_start >= end_idx) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, end_idx);

		for (uint32_t vidx = batch_start; vidx < batch_end; vidx++) {
			uint32_t missing_ct = 0;

			if (gstate.need_missingness && lstate.initialized) {
				auto *missingness = lstate.missingness_buf.As<uintptr_t>();
				auto *genovec = lstate.genovec_buf.As<uintptr_t>();

				plink2::PglErr err = plink2::PgrGetMissingness(sample_include, lstate.pssi, sample_ct, vidx,
				                                               &lstate.pgr, missingness, genovec);

				if (err != plink2::kPglRetSuccess) {
					throw IOException("plink_missing: PgrGetMissingness failed for variant %u", vidx);
				}

				missing_ct = static_cast<uint32_t>(plink2::PopcountWords(missingness, word_ct));
			}

			uint32_t obs_ct = sample_ct - missing_ct;
			double f_miss = (sample_ct > 0) ? static_cast<double>(missing_ct) / static_cast<double>(sample_ct) : 0.0;

			// Fill projected columns
			for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
				auto file_col = column_ids[out_col];
				if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
					continue;
				}

				auto &vec = output.data[out_col];

				switch (file_col) {
				case VCOL_CHROM: {
					auto &val = bind_data.variants.chroms[vidx];
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case VCOL_POS: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.positions[vidx];
					break;
				}
				case VCOL_ID: {
					auto &val = bind_data.variants.ids[vidx];
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case VCOL_REF: {
					auto &val = bind_data.variants.refs[vidx];
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case VCOL_ALT: {
					auto &val = bind_data.variants.alts[vidx];
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case VCOL_MISSING_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(missing_ct);
					break;
				}
				case VCOL_OBS_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(obs_ct);
					break;
				}
				case VCOL_F_MISS: {
					FlatVector::GetData<double>(vec)[rows_emitted] = f_miss;
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
// Scan — sample mode (single-threaded: scan all variants, then emit rows)
// ---------------------------------------------------------------------------

static void PlinkMissingScanSample(const PlinkMissingBindData &bind_data, PlinkMissingGlobalState &gstate,
                                   PlinkMissingLocalState &lstate, DataChunk &output) {
	uint32_t sample_ct = bind_data.effective_sample_ct;

	// Phase 1: Scan all variants, accumulate per-sample missing counts
	if (!gstate.variant_scan_done && gstate.need_missingness && lstate.initialized) {
		uint32_t end_idx = gstate.end_variant_idx;

		const uintptr_t *sample_include = nullptr;
		if (bind_data.has_sample_subset && bind_data.sample_subset) {
			sample_include = bind_data.sample_subset->SampleInclude();
		}

		auto *missingness = lstate.missingness_buf.As<uintptr_t>();
		auto *genovec = lstate.genovec_buf.As<uintptr_t>();
		uintptr_t word_ct = plink2::BitCtToWordCt(sample_ct);

		for (uint32_t vidx = gstate.start_variant_idx; vidx < end_idx; vidx++) {
			plink2::PglErr err = plink2::PgrGetMissingness(sample_include, lstate.pssi, sample_ct, vidx, &lstate.pgr,
			                                               missingness, genovec);

			if (err != plink2::kPglRetSuccess) {
				throw IOException("plink_missing: PgrGetMissingness failed for variant %u", vidx);
			}

			// Iterate set bits (missing samples) and increment per-sample counters
			for (uintptr_t w = 0; w < word_ct; w++) {
				uintptr_t word = missingness[w];
				while (word) {
					uint32_t bit_pos = plink2::ctzw(word);
					uint32_t sample_idx = static_cast<uint32_t>(w) * plink2::kBitsPerWord + bit_pos;
					if (sample_idx < sample_ct) {
						gstate.sample_missing_counts[sample_idx]++;
					}
					word &= word - 1; // clear lowest set bit
				}
			}
		}

		gstate.variant_scan_done = true;
	}

	// Phase 2: Emit sample rows from accumulated counts
	auto &column_ids = gstate.column_ids;
	uint32_t total_variant_ct = gstate.total_variant_ct;

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t sidx = gstate.next_sample_idx.fetch_add(1);
		if (sidx >= sample_ct) {
			break;
		}

		uint32_t missing_ct = gstate.need_missingness ? gstate.sample_missing_counts[sidx] : 0;
		uint32_t obs_ct = total_variant_ct - missing_ct;
		double f_miss =
		    (total_variant_ct > 0) ? static_cast<double>(missing_ct) / static_cast<double>(total_variant_ct) : 0.0;

		// Map subset position → original sample index for FID/IID lookup
		uint32_t orig_idx = bind_data.has_sample_subset ? bind_data.subset_original_indices[sidx] : sidx;

		for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
			auto file_col = column_ids[out_col];
			if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
				continue;
			}

			auto &vec = output.data[out_col];

			switch (file_col) {
			case SCOL_FID: {
				if (bind_data.has_sample_info && !bind_data.sample_info.fids.empty() &&
				    orig_idx < bind_data.sample_info.fids.size() && !bind_data.sample_info.fids[orig_idx].empty()) {
					FlatVector::GetData<string_t>(vec)[rows_emitted] =
					    StringVector::AddString(vec, bind_data.sample_info.fids[orig_idx]);
				} else {
					FlatVector::SetNull(vec, rows_emitted, true);
				}
				break;
			}
			case SCOL_IID: {
				FlatVector::GetData<string_t>(vec)[rows_emitted] =
				    StringVector::AddString(vec, bind_data.sample_info.iids[orig_idx]);
				break;
			}
			case SCOL_MISSING_CT: {
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(missing_ct);
				break;
			}
			case SCOL_OBS_CT: {
				FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(obs_ct);
				break;
			}
			case SCOL_F_MISS: {
				FlatVector::GetData<double>(vec)[rows_emitted] = f_miss;
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
// Scan dispatcher
// ---------------------------------------------------------------------------

static void PlinkMissingScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkMissingBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkMissingGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkMissingLocalState>();

	if (gstate.sample_mode) {
		PlinkMissingScanSample(bind_data, gstate, lstate, output);
	} else {
		PlinkMissingScanVariant(bind_data, gstate, lstate, output);
	}
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkMissing(ExtensionLoader &loader) {
	TableFunction plink_missing("plink_missing", {LogicalType::VARCHAR}, PlinkMissingScan, PlinkMissingBind,
	                            PlinkMissingInitGlobal, PlinkMissingInitLocal);

	plink_missing.projection_pushdown = true;

	plink_missing.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_missing.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_missing.named_parameters["mode"] = LogicalType::VARCHAR;
	plink_missing.named_parameters["samples"] = LogicalType::ANY;
	plink_missing.named_parameters["region"] = LogicalType::VARCHAR;

	loader.RegisterFunction(plink_missing);
}

} // namespace duckdb
