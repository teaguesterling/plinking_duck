#include "plink_freq.hpp"
#include "plink_common.hpp"

#include <atomic>

namespace duckdb {

// ---------------------------------------------------------------------------
// Column indices
// ---------------------------------------------------------------------------

// Default columns: CHROM(0) POS(1) ID(2) REF(3) ALT(4) ALT_FREQ(5) OBS_CT(6)
// With counts := true: + HOM_REF_CT(7) HET_CT(8) HOM_ALT_CT(9) MISSING_CT(10)
static constexpr idx_t COL_CHROM = 0;
static constexpr idx_t COL_POS = 1;
static constexpr idx_t COL_ID = 2;
static constexpr idx_t COL_REF = 3;
static constexpr idx_t COL_ALT = 4;
static constexpr idx_t COL_ALT_FREQ = 5;
static constexpr idx_t COL_OBS_CT = 6;
static constexpr idx_t COL_HOM_REF_CT = 7;
static constexpr idx_t COL_HET_CT = 8;
static constexpr idx_t COL_HOM_ALT_CT = 9;
static constexpr idx_t COL_MISSING_CT = 10;

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkFreqBindData : public TableFunctionData {
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
	uint32_t effective_sample_ct = 0; // subset_sample_ct or raw_sample_ct

	// Region filtering
	VariantRange variant_range;

	// Options
	bool include_counts = false;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkFreqGlobalState : public GlobalTableFunctionState {
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t start_variant_idx = 0;
	uint32_t end_variant_idx = 0;
	vector<column_t> column_ids;
	bool need_frequencies = false; // true if any freq/count column is projected

	idx_t MaxThreads() const override {
		uint32_t range = end_variant_idx - start_variant_idx;
		return std::min<idx_t>(range / 500 + 1, 16);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkFreqLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	bool initialized = false;

	~PlinkFreqLocalState() {
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

static unique_ptr<FunctionData> PlinkFreqBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkFreqBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "counts") {
			bind_data->include_counts = kv.second.GetValue<bool>();
		} else if (kv.first == "dosage") {
			if (kv.second.GetValue<bool>()) {
				throw NotImplementedException("plink_freq: dosage mode is not yet implemented");
			}
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled after pgenlib init
		}
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_freq: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		// .psam is optional for plink_freq — frequency computation only needs sample count
	}

	// --- Initialize pgenlib (Phase 1) to get counts ---
	plink2::PgenFileInfo pgfi;
	plink2::PreinitPgfi(&pgfi);

	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr,
	                                            UINT32_MAX, UINT32_MAX,
	                                            &header_ctrl, &pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("plink_freq: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;

	// Phase 2 (needed to validate the file, but per-thread readers will re-init)
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
		throw IOException("plink_freq: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "plink_freq");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_freq: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (optional) ---
	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);
		bind_data->has_sample_info = true;

		if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_freq: sample count mismatch: .pgen has %u samples, "
			                            ".psam/.fam '%s' has %llu samples",
			                            bind_data->raw_sample_ct, bind_data->psam_path,
			                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
		}
	}

	// --- Process samples parameter ---
	bind_data->effective_sample_ct = bind_data->raw_sample_ct;

	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		auto indices = ResolveSampleIndices(
		    samples_it->second, bind_data->raw_sample_ct,
		    bind_data->has_sample_info ? &bind_data->sample_info : nullptr, "plink_freq");

		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range = ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_freq");
	}

	// --- Register output columns ---
	names = {"CHROM", "POS", "ID", "REF", "ALT", "ALT_FREQ", "OBS_CT"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::INTEGER};

	if (bind_data->include_counts) {
		names.insert(names.end(), {"HOM_REF_CT", "HET_CT", "HOM_ALT_CT", "MISSING_CT"});
		return_types.insert(return_types.end(),
		                    {LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER});
	}

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkFreqInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkFreqBindData>();
	auto state = make_uniq<PlinkFreqGlobalState>();

	if (bind_data.variant_range.has_filter) {
		state->start_variant_idx = bind_data.variant_range.start_idx;
		state->end_variant_idx = bind_data.variant_range.end_idx;
	} else {
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.raw_variant_ct;
	}

	state->next_variant_idx.store(state->start_variant_idx);
	state->column_ids = input.column_ids;

	// Check if any frequency/count columns are projected
	state->need_frequencies = false;
	for (auto col_id : input.column_ids) {
		if (col_id != COLUMN_IDENTIFIER_ROW_ID && col_id >= COL_ALT_FREQ && col_id <= COL_MISSING_CT) {
			state->need_frequencies = true;
			break;
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkFreqInitLocal(ExecutionContext &context,
                                                               TableFunctionInitInput &input,
                                                               GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkFreqBindData>();
	auto &gstate = global_state->Cast<PlinkFreqGlobalState>();
	auto state = make_uniq<PlinkFreqLocalState>();

	if (!gstate.need_frequencies) {
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
		throw IOException("plink_freq: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_freq: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_freq: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static constexpr uint32_t FREQ_BATCH_SIZE = 128;

static void PlinkFreqScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkFreqBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkFreqGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkFreqLocalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t end_idx = gstate.end_variant_idx;
	uint32_t sample_ct = bind_data.effective_sample_ct;

	// Get sample subsetting pointers (nullptr if no subset)
	const uintptr_t *sample_include = nullptr;
	const uintptr_t *interleaved_vec = nullptr;
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		sample_include = bind_data.sample_subset->SampleInclude();
		interleaved_vec = bind_data.sample_subset->InterleavedVec();
	}

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(FREQ_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
		if (batch_start >= end_idx) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, end_idx);

		for (uint32_t vidx = batch_start; vidx < batch_end; vidx++) {

			// Compute genotype counts using PgrGetCounts fast-path
			STD_ARRAY_DECL(uint32_t, 4, genocounts);
			genocounts[0] = genocounts[1] = genocounts[2] = genocounts[3] = 0;

			if (gstate.need_frequencies && lstate.initialized) {
				plink2::PglErr err = plink2::PgrGetCounts(sample_include, interleaved_vec, lstate.pssi, sample_ct,
				                                          vidx, &lstate.pgr, genocounts);

				if (err != plink2::kPglRetSuccess) {
					throw IOException("plink_freq: PgrGetCounts failed for variant %u", vidx);
				}
			}

			// Compute derived values
			uint32_t obs_sample_ct = genocounts[0] + genocounts[1] + genocounts[2];
			uint32_t obs_ct = 2 * obs_sample_ct; // allele observations
			double alt_freq;
			bool freq_is_null = false;

			if (obs_sample_ct == 0) {
				freq_is_null = true;
				alt_freq = 0.0;
			} else {
				alt_freq = (static_cast<double>(genocounts[1]) + 2.0 * static_cast<double>(genocounts[2])) /
				           (2.0 * static_cast<double>(obs_sample_ct));
			}

			// Fill projected columns
			for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
				auto file_col = column_ids[out_col];
				if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
					continue;
				}

				auto &vec = output.data[out_col];

				switch (file_col) {
				case COL_CHROM: {
					auto &val = bind_data.variants.chroms[vidx];
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case COL_POS: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.positions[vidx];
					break;
				}
				case COL_ID: {
					auto &val = bind_data.variants.ids[vidx];
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_REF: {
					auto &val = bind_data.variants.refs[vidx];
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case COL_ALT: {
					auto &val = bind_data.variants.alts[vidx];
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_ALT_FREQ: {
					if (freq_is_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = alt_freq;
					}
					break;
				}
				case COL_OBS_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(obs_ct);
					break;
				}
				case COL_HOM_REF_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(genocounts[0]);
					break;
				}
				case COL_HET_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(genocounts[1]);
					break;
				}
				case COL_HOM_ALT_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(genocounts[2]);
					break;
				}
				case COL_MISSING_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(genocounts[3]);
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
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkFreq(ExtensionLoader &loader) {
	TableFunction plink_freq("plink_freq", {LogicalType::VARCHAR}, PlinkFreqScan, PlinkFreqBind, PlinkFreqInitGlobal,
	                         PlinkFreqInitLocal);

	plink_freq.projection_pushdown = true;

	plink_freq.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_freq.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_freq.named_parameters["samples"] = LogicalType::ANY;
	plink_freq.named_parameters["region"] = LogicalType::VARCHAR;
	plink_freq.named_parameters["counts"] = LogicalType::BOOLEAN;
	plink_freq.named_parameters["dosage"] = LogicalType::BOOLEAN;

	loader.RegisterFunction(plink_freq);
}

} // namespace duckdb
