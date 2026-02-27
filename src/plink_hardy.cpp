#include "plink_hardy.hpp"
#include "plink_common.hpp"

#include <atomic>
#include <cmath>

namespace duckdb {

// ---------------------------------------------------------------------------
// Column indices
// ---------------------------------------------------------------------------

// CHROM(0) POS(1) ID(2) REF(3) ALT(4) A1(5) HOM_REF_CT(6) HET_CT(7)
// HOM_ALT_CT(8) O_HET(9) E_HET(10) P_HWE(11)
static constexpr idx_t COL_CHROM = 0;
static constexpr idx_t COL_POS = 1;
static constexpr idx_t COL_ID = 2;
static constexpr idx_t COL_REF = 3;
static constexpr idx_t COL_ALT = 4;
static constexpr idx_t COL_A1 = 5;
static constexpr idx_t COL_HOM_REF_CT = 6;
static constexpr idx_t COL_HET_CT = 7;
static constexpr idx_t COL_HOM_ALT_CT = 8;
static constexpr idx_t COL_O_HET = 9;
static constexpr idx_t COL_E_HET = 10;
static constexpr idx_t COL_P_HWE = 11;

// ---------------------------------------------------------------------------
// HWE exact test (Wigginton et al. 2005)
// ---------------------------------------------------------------------------

// Computes the Hardy-Weinberg equilibrium exact test p-value.
// Algorithm: enumerate all valid heterozygote counts for fixed allele counts,
// compute relative probabilities via a recurrence relation, then sum
// probabilities of configurations as likely or less likely than observed.
static double HweExactTest(int32_t obs_hom1, int32_t obs_hets, int32_t obs_hom2, bool midp) {
	if (obs_hom1 + obs_hets + obs_hom2 == 0) {
		return 1.0;
	}

	// Order so obs_homr <= obs_homc (rare/common homozygotes)
	int32_t obs_homc = std::max(obs_hom1, obs_hom2);
	int32_t obs_homr = std::min(obs_hom1, obs_hom2);

	int32_t rare_copies = 2 * obs_homr + obs_hets;
	int32_t common_copies = 2 * obs_homc + obs_hets;
	int32_t n = obs_homc + obs_homr + obs_hets;

	// het counts must have same parity as rare_copies
	// Find the mode: expected het count under HWE
	int32_t mid = static_cast<int32_t>(static_cast<double>(rare_copies) * common_copies / (2.0 * n));
	// Ensure same parity as rare_copies
	if ((mid % 2) != (rare_copies % 2)) {
		mid++;
	}

	// Allocate probability array for valid het counts (0..rare_copies, stepping by 2)
	// We store all entries but only use those with correct parity
	vector<double> het_probs(rare_copies + 1, 0.0);
	het_probs[mid] = 1.0;
	double sum = 1.0;

	// Fill upward from mid (increasing het count)
	int32_t curr_hets = mid;
	int32_t curr_homr = (rare_copies - mid) / 2;
	int32_t curr_homc = (common_copies - mid) / 2;

	while (curr_hets <= rare_copies - 2) {
		// P(k+2)/P(k) = 4 * homr * homc / ((k+1) * (k+2))
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /
		                           ((static_cast<double>(curr_hets) + 1.0) * (static_cast<double>(curr_hets) + 2.0));
		sum += het_probs[curr_hets + 2];
		curr_homr--;
		curr_homc--;
		curr_hets += 2;
	}

	// Fill downward from mid (decreasing het count)
	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = (common_copies - mid) / 2;

	while (curr_hets >= 2) {
		// P(k-2)/P(k) = k * (k-1) / (4 * (homr+1) * (homc+1))
		het_probs[curr_hets - 2] =
		    het_probs[curr_hets] * static_cast<double>(curr_hets) * (static_cast<double>(curr_hets) - 1.0) /
		    (4.0 * (static_cast<double>(curr_homr) + 1.0) * (static_cast<double>(curr_homc) + 1.0));
		sum += het_probs[curr_hets - 2];
		curr_homr++;
		curr_homc++;
		curr_hets -= 2;
	}

	// Sum probabilities of all het counts with P <= P(observed)
	double obs_prob = het_probs[obs_hets] / sum;
	double p_value = 0.0;

	// Use a small tolerance for floating-point comparison
	double threshold = obs_prob * (1.0 + 1e-8);

	for (int32_t i = 0; i <= rare_copies; i += 2) {
		if (het_probs[i] / sum <= threshold) {
			p_value += het_probs[i] / sum;
		}
	}
	// Also check odd values if rare_copies parity is odd
	if (rare_copies % 2 == 1) {
		for (int32_t i = 1; i <= rare_copies; i += 2) {
			if (het_probs[i] / sum <= threshold) {
				p_value += het_probs[i] / sum;
			}
		}
	}

	if (midp) {
		p_value -= obs_prob * 0.5;
	}

	if (p_value < 0.0) {
		return 0.0;
	}
	if (p_value > 1.0) {
		return 1.0;
	}
	return p_value;
}

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkHardyBindData : public TableFunctionData {
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

	// Region filtering
	VariantRange variant_range;

	// Options
	bool midp = false;
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkHardyGlobalState : public GlobalTableFunctionState {
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t start_variant_idx = 0;
	uint32_t end_variant_idx = 0;
	vector<column_t> column_ids;
	bool need_genotype_counts = false;

	idx_t MaxThreads() const override {
		uint32_t range = end_variant_idx - start_variant_idx;
		return std::min<idx_t>(range / 500 + 1, 16);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkHardyLocalState : public LocalTableFunctionState {
	plink2::PgenFileInfo pgfi;
	AlignedBuffer pgfi_alloc_buf;

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	bool initialized = false;

	~PlinkHardyLocalState() {
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

static unique_ptr<FunctionData> PlinkHardyBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkHardyBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "midp") {
			bind_data->midp = kv.second.GetValue<bool>();
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled after pgenlib init
		}
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_hardy: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		// .psam is optional for plink_hardy — HWE only needs sample count
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
		throw IOException("plink_hardy: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;

	// Phase 2 (validates the file; per-thread readers will re-init)
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
		throw IOException("plink_hardy: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "plink_hardy");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_hardy: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (optional) ---
	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);
		bind_data->has_sample_info = true;

		if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_hardy: sample count mismatch: .pgen has %u samples, "
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
		                         bind_data->has_sample_info ? &bind_data->sample_info : nullptr, "plink_hardy");

		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range =
		    ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_hardy");
	}

	// --- Register output columns ---
	names = {"CHROM", "POS", "ID", "REF", "ALT", "A1", "HOM_REF_CT", "HET_CT", "HOM_ALT_CT", "O_HET", "E_HET", "P_HWE"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	                LogicalType::INTEGER, LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE};

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkHardyInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkHardyBindData>();
	auto state = make_uniq<PlinkHardyGlobalState>();

	if (bind_data.variant_range.has_filter) {
		state->start_variant_idx = bind_data.variant_range.start_idx;
		state->end_variant_idx = bind_data.variant_range.end_idx;
	} else {
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.raw_variant_ct;
	}

	state->next_variant_idx.store(state->start_variant_idx);
	state->column_ids = input.column_ids;

	// Check if any genotype-dependent columns are projected
	state->need_genotype_counts = false;
	for (auto col_id : input.column_ids) {
		if (col_id != COLUMN_IDENTIFIER_ROW_ID && col_id >= COL_HOM_REF_CT && col_id <= COL_P_HWE) {
			state->need_genotype_counts = true;
			break;
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkHardyInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                               GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkHardyBindData>();
	auto &gstate = global_state->Cast<PlinkHardyGlobalState>();
	auto state = make_uniq<PlinkHardyLocalState>();

	if (!gstate.need_genotype_counts) {
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
		throw IOException("plink_hardy: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_hardy: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_hardy: PgrInit failed for '%s'", bind_data.pgen_path);
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

static constexpr uint32_t HARDY_BATCH_SIZE = 128;

static void PlinkHardyScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkHardyBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkHardyGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkHardyLocalState>();

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
		uint32_t claim_size = std::min(HARDY_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
		if (batch_start >= end_idx) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, end_idx);

		for (uint32_t vidx = batch_start; vidx < batch_end; vidx++) {

			// Compute genotype counts using PgrGetCounts fast-path
			STD_ARRAY_DECL(uint32_t, 4, genocounts);
			genocounts[0] = genocounts[1] = genocounts[2] = genocounts[3] = 0;

			if (gstate.need_genotype_counts && lstate.initialized) {
				plink2::PglErr err = plink2::PgrGetCounts(sample_include, interleaved_vec, lstate.pssi, sample_ct, vidx,
				                                          &lstate.pgr, genocounts);

				if (err != plink2::kPglRetSuccess) {
					throw IOException("plink_hardy: PgrGetCounts failed for variant %u", vidx);
				}
			}

			// Compute derived values
			uint32_t hom_ref = genocounts[0];
			uint32_t het = genocounts[1];
			uint32_t hom_alt = genocounts[2];
			uint32_t obs = hom_ref + het + hom_alt;

			double o_het;
			double e_het;
			double p_hwe;
			bool stats_are_null = (obs == 0);

			if (stats_are_null) {
				o_het = 0.0;
				e_het = 0.0;
				p_hwe = 1.0;
			} else {
				o_het = static_cast<double>(het) / static_cast<double>(obs);
				double p = (2.0 * hom_ref + het) / (2.0 * obs);
				double q = 1.0 - p;
				e_het = 2.0 * p * q;
				p_hwe = HweExactTest(static_cast<int32_t>(hom_ref), static_cast<int32_t>(het),
				                     static_cast<int32_t>(hom_alt), bind_data.midp);
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
				case COL_A1: {
					// A1 = tested allele (alternate)
					auto &val = bind_data.variants.alts[vidx];
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_HOM_REF_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(hom_ref);
					break;
				}
				case COL_HET_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(het);
					break;
				}
				case COL_HOM_ALT_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(hom_alt);
					break;
				}
				case COL_O_HET: {
					if (stats_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = o_het;
					}
					break;
				}
				case COL_E_HET: {
					if (stats_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = e_het;
					}
					break;
				}
				case COL_P_HWE: {
					if (stats_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = p_hwe;
					}
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

void RegisterPlinkHardy(ExtensionLoader &loader) {
	TableFunction plink_hardy("plink_hardy", {LogicalType::VARCHAR}, PlinkHardyScan, PlinkHardyBind,
	                          PlinkHardyInitGlobal, PlinkHardyInitLocal);

	plink_hardy.projection_pushdown = true;

	plink_hardy.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_hardy.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_hardy.named_parameters["samples"] = LogicalType::ANY;
	plink_hardy.named_parameters["region"] = LogicalType::VARCHAR;
	plink_hardy.named_parameters["midp"] = LogicalType::BOOLEAN;

	loader.RegisterFunction(plink_hardy);
}

} // namespace duckdb
