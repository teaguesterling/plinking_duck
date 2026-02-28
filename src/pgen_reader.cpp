#include "pgen_reader.hpp"
#include "plink_common.hpp"

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

	// Sample subsetting
	bool has_sample_subset = false;
	vector<uint32_t> sample_indices; // 0-based indices into .pgen sample order
	uint32_t subset_sample_ct = 0;
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
		} else if (kv.first == "samples") {
			// Handled after pgenlib init (need raw_sample_ct)
		}
	}

	if (bind_data->include_dosages) {
		throw NotImplementedException("read_pgen: dosages support is not yet implemented");
	}
	if (bind_data->include_phased) {
		throw NotImplementedException("read_pgen: phased support is not yet implemented");
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

	// --- Register output columns ---
	names = {"CHROM", "POS", "ID", "REF", "ALT", "genotypes"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::TINYINT)};

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

	// Check if genotypes column is in the projection
	state->need_genotypes = false;
	for (auto col_id : input.column_ids) {
		if (col_id == PgenBindData::GENOTYPES_COL_IDX) {
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

				plink2::PglErr err = plink2::PgrGet(sample_include, lstate.pssi, output_sample_ct, vidx, &lstate.pgr,
				                                    lstate.genovec_buf.As<uintptr_t>());

				if (err != plink2::kPglRetSuccess) {
					throw IOException("read_pgen: PgrGet failed for variant %u", vidx);
				}

				// Unpack 2-bit genotypes to bytes (-9 = missing)
				plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), output_sample_ct,
				                             lstate.genotype_bytes.data());
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
				case 5: { // genotypes — LIST(TINYINT)
					if (!genotypes_read) {
						FlatVector::SetNull(vec, rows_emitted, true);
						break;
					}

					// Build the LIST(TINYINT) entry
					auto list_size = static_cast<idx_t>(output_sample_ct);
					auto current_offset = ListVector::GetListSize(vec);
					auto &entry = FlatVector::GetData<list_entry_t>(vec)[rows_emitted];
					entry.offset = current_offset;
					entry.length = list_size;

					ListVector::Reserve(vec, current_offset + list_size);
					auto &child = ListVector::GetEntry(vec);
					auto *child_data = FlatVector::GetData<int8_t>(child);
					auto &child_validity = FlatVector::Validity(child);

					for (idx_t s = 0; s < list_size; s++) {
						int8_t geno = lstate.genotype_bytes[s];
						if (geno == -9) {
							child_validity.SetInvalid(current_offset + s);
							child_data[current_offset + s] = 0; // placeholder
						} else {
							child_data[current_offset + s] = geno;
						}
					}

					ListVector::SetListSize(vec, current_offset + list_size);
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

	loader.RegisterFunction(read_pgen);
}

} // namespace duckdb
