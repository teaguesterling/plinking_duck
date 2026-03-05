#include "pgen_reader.hpp"
#include "plink_common.hpp"

#include "duckdb/common/string_util.hpp"

#include <algorithm>
#include <atomic>

namespace duckdb {

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PgenBindData : public TableFunctionData {
	string pgen_path;
	string pvar_path;
	string psam_path;

	// Variant metadata (offset-indexed for on-demand parsing)
	VariantMetadataIndex variants;

	// Sample metadata (optional — may not have .psam)
	SampleInfo sample_info;
	bool has_sample_info = false;
	uint32_t sample_ct = 0; // from .pgen header if no .psam

	// pgenlib initialization results (needed by init/scan)
	uint32_t raw_variant_ct = 0;
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

	// Columns mode layout (genotypes := 'columns')
	vector<string> genotype_column_names;     // IIDs for column names
	idx_t columns_mode_first_geno_col = 0;    // first genotype column index in schema
	uint32_t columns_mode_geno_col_count = 0; // number of genotype columns
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PgenGlobalState : public GlobalTableFunctionState {
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t total_variants = 0;

	// Projection flags (computed once in init)
	bool need_genotypes = false;
	vector<column_t> column_ids;

	idx_t MaxThreads() const override {
		return std::min<idx_t>(total_variants / 1000 + 1, 16);
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
// Bind function
// ---------------------------------------------------------------------------

static unique_ptr<FunctionData> PgenBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PgenBindData>();
	bind_data->pgen_path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
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
		}
	}

	if (bind_data->include_dosages && bind_data->include_phased) {
		throw InvalidInputException("read_pgen: dosages and phased cannot both be true");
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFile(fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("read_pgen: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFile(fs, bind_data->pgen_path, {".psam", ".fam"});
		// .psam is optional for read_pgen — if not found, we operate in index-only mode
	}

	// --- Initialize pgenlib (Phase 1) ---
	plink2::PgenFileInfo pgfi;
	plink2::PreinitPgfi(&pgfi);

	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr, // no .pgi file
	                                            UINT32_MAX,                            // infer raw_variant_ct
	                                            UINT32_MAX,                            // infer raw_sample_ct
	                                            &header_ctrl, &pgfi, &pgfi_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&pgfi, &cleanup_err);
		throw IOException("read_pgen: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;

	// --- Phase 2: allocate pgfi memory and complete init ---
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
		throw IOException("read_pgen: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadataIndex(context, bind_data->pvar_path, "read_pgen");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("read_pgen: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (optional) ---
	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleInfo(context, bind_data->psam_path);
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

	// --- Resolve genotype output mode ---
	uint32_t output_sample_ct = bind_data->has_sample_subset ? bind_data->subset_sample_ct : bind_data->sample_ct;

	string genotypes_str = "auto";
	auto genotypes_it = input.named_parameters.find("genotypes");
	if (genotypes_it != input.named_parameters.end()) {
		genotypes_str = genotypes_it->second.GetValue<string>();
	}
	bind_data->genotype_mode = ResolveGenotypeMode(genotypes_str, output_sample_ct, "read_pgen");

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

static unique_ptr<GlobalTableFunctionState> PgenInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PgenBindData>();
	auto state = make_uniq<PgenGlobalState>();

	state->total_variants = bind_data.raw_variant_ct;
	state->column_ids = input.column_ids;

	// Check if genotypes column(s) are in the projection
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
		} else if (col_id == PgenBindData::GENOTYPES_COL_IDX) {
			state->need_genotypes = true;
			break;
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

	if (!gstate.need_genotypes) {
		// No genotype columns needed — skip pgenlib initialization entirely
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
		throw IOException("read_pgen: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("read_pgen: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("read_pgen: PgrInit failed for '%s'", bind_data.pgen_path);
	}

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

	// Set up sample subsetting if needed
	if (bind_data.has_sample_subset) {
		// Build sample_include bitmask
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

		plink2::PgrSetSampleSubsetIndex(cumulative_popcounts, &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	state->initialized = true;
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

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		// Claim a batch of variants, capped to remaining output capacity
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(PGEN_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
		if (batch_start >= total_variants) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, total_variants);

		for (uint32_t vidx = batch_start; vidx < batch_end; vidx++) {

			// Read genotype data if needed (before filling columns, since
			// we need it for the genotypes column)
			bool genotypes_read = false;
			if (gstate.need_genotypes && lstate.initialized) {
				const uintptr_t *sample_include =
				    bind_data.has_sample_subset ? lstate.sample_include_buf.As<uintptr_t>() : nullptr;

				if (bind_data.include_dosages) {
					uint32_t dosage_ct = 0;
					plink2::PglErr err = plink2::PgrGetD(
					    sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
					    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
					    lstate.dosage_main_buf.As<uint16_t>(), &dosage_ct);
					if (err != plink2::kPglRetSuccess) {
						throw IOException("read_pgen: PgrGetD failed for variant %u", vidx);
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
						throw IOException("read_pgen: PgrGetP failed for variant %u", vidx);
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
					auto val = bind_data.variants.GetChrom(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case 1: { // POS
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.GetPos(vidx);
					break;
				}
				case 2: { // ID
					auto val = bind_data.variants.GetId(vidx);
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case 3: { // REF
					auto val = bind_data.variants.GetRef(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case 4: { // ALT
					auto val = bind_data.variants.GetAlt(vidx);
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
								if (geno == -9) {
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
								if (geno == -9) {
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
							if (bind_data.include_dosages) {
								double dosage = lstate.dosage_doubles[sample_pos];
								if (dosage == -9.0) {
									FlatVector::SetNull(vec, rows_emitted, true);
								} else {
									FlatVector::GetData<double>(vec)[rows_emitted] = dosage;
								}
							} else {
								int8_t geno = lstate.genotype_bytes[sample_pos];
								if (geno == -9) {
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
// Registration
// ---------------------------------------------------------------------------

void RegisterPgenReader(ExtensionLoader &loader) {
	TableFunction read_pgen("read_pgen", {LogicalType::VARCHAR}, PgenScan, PgenBind, PgenInitGlobal, PgenInitLocal);

	read_pgen.projection_pushdown = true;

	read_pgen.named_parameters["pvar"] = LogicalType::VARCHAR;
	read_pgen.named_parameters["psam"] = LogicalType::VARCHAR;
	read_pgen.named_parameters["dosages"] = LogicalType::BOOLEAN;
	read_pgen.named_parameters["phased"] = LogicalType::BOOLEAN;
	// Accept ANY for samples — type dispatch (LIST(INTEGER) vs LIST(VARCHAR))
	// is handled in PgenBind based on the actual value type.
	read_pgen.named_parameters["samples"] = LogicalType::ANY;
	read_pgen.named_parameters["genotypes"] = LogicalType::VARCHAR;
	read_pgen.named_parameters["orient"] = LogicalType::VARCHAR;

	loader.RegisterFunction(read_pgen);
}

} // namespace duckdb
