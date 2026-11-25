#include "plink_freq.hpp"
#include "pgen_header_cache.hpp"
#include "duckdb_compat.hpp"
#include "register_helper.hpp"
#include "named_parameter_compat.hpp" // DeclareNamedParameters / NamedParamMap (DuckDB v2.0 typed kwargs)
#include "plink_common.hpp"
#include "pgen_vfs_opener.hpp"
#include "catalog_cohort.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>

namespace duckdb {

// ---------------------------------------------------------------------------
// Column indices
// ---------------------------------------------------------------------------

// Default columns: CHROM(0) POS(1) ID(2) REF(3) ALT(4) ALT_FREQ(5) OBS_CT(6)
// With counts := true: + HOM_REF_CT(7) HET_CT(8) HOM_ALT_CT(9) MISSING_CT(10)
// With dosage := true: + IMP_R2(11) (only present when dosage is enabled)
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
static constexpr idx_t COL_IMP_R2 = 11;

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkFreqBindData : public TableFunctionData {
	string pgen_path;
	bool use_vfs = false;             // route .pgen opens through DuckDB's VFS (plinking_pgen_io)
	PgenLocalizeGuard localize_guard; // owns downloaded temp .pgen for 'localize' (per-query)
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
	uint32_t effective_sample_ct = 0; // subset_sample_ct or raw_sample_ct

	// Region filtering
	VariantRange variant_range;

	// Options
	bool include_counts = false;
	bool include_dosage = false;
	bool file_has_dosage = false; // true if pgen file contains dosage data

	// Ploidy/sex-aware handling for chrX/Y/MT
	ParBounds par_bounds;        // pseudo-autosomal boundaries for the genome build
	vector<uint8_t> aligned_sex; // sex per sample in effective (post-subset) order
	bool have_sex = false;       // true iff aligned_sex is populated

	// Dynamic column index for IMP_R2 (depends on whether counts is enabled)
	idx_t imp_r2_col_idx = 0;

	// A catalog read (cohort set): the scan reads the cohort's segments through a ConcatPgenReader
	// per variant block. raw_sample_ct is then the cohort's sample count, and the subset is held as
	// sorted cohort-wide indices (ConcatPgenReader splits it per sample block).
	unique_ptr<CatalogCohort> cohort;
	bool cohort_has_subset = false;
	vector<uint32_t> cohort_subset_sorted;
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
	uint32_t max_threads_config = 0;

	idx_t MaxThreads() const override {
		uint32_t range = end_variant_idx - start_variant_idx;
		return ApplyMaxThreadsCap(range / 500 + 1, max_threads_config);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkFreqLocalState : public LocalTableFunctionState {
	shared_ptr<PgenHeader> header; //! decoded .pgen header (cache tier 1); pgr points into it

	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	// Raw-hardcall buffers for the ploidy/sex-aware path (chrX/Y/MT only).
	AlignedBuffer genovec_buf;
	vector<int8_t> geno_bytes;

	bool initialized = false;

	// Catalog read: the reader of the variant block last read.
	CatalogCohortReader cohort_reader;
	bool cohort_ready = false;

	~PlinkFreqLocalState() override {
		if (initialized) {
			plink2::PglErr reterr = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr, &reterr);
			// the header is released after pgr (member order)
		}
	}
};

// ---------------------------------------------------------------------------
// Bind
// ---------------------------------------------------------------------------

static void SetFreqColumns(PlinkFreqBindData &bind_data, vector<LogicalType> &return_types, vector<CompatName> &names);

//! Bind over a pfile catalog: every variant block in catalog order, counted over every sample block.
//! samples, region, selection, counts and build work as over one fileset. dosage := true does not:
//! ConcatPgenReader has no dosage count (pvar and psam are rejected by ResolveAnalysisCatalog).
static void BindFreqCatalog(ClientContext &context, TableFunctionBindInput &input, PlinkFreqBindData &bind_data,
                            const string &catalog_path) {
	if (bind_data.include_dosage) {
		throw InvalidInputException("plink_freq: dosage := true is not supported with a catalog (%s)", catalog_path);
	}
	bind_data.cohort = LoadCatalogCohort(context, catalog_path, "plink_freq");
	auto &cohort = *bind_data.cohort;
	Selection selection = ParseSelection(input.named_parameters);
	ResolvedSelection resolved = ResolveCohortSelection(selection, cohort, "plink_freq");
	bind_data.raw_sample_ct = cohort.sample_ct;
	bind_data.raw_variant_ct = cohort.variant_ct;
	bind_data.effective_sample_ct = resolved.effective_sample_ct;
	bind_data.cohort_has_subset = resolved.has_subset;
	bind_data.cohort_subset_sorted = resolved.subset_sorted;
	bind_data.has_sample_info = true;
	if (!cohort.samples.sexes.empty()) {
		bind_data.aligned_sex =
		    BuildAlignedSex(cohort.samples, resolved.has_subset ? &resolved.subset_sorted : nullptr);
		bind_data.have_sex = !bind_data.aligned_sex.empty();
	}
}

static unique_ptr<FunctionData> PlinkFreqBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<CompatName> &names) {
	auto bind_data = make_uniq<PlinkFreqBindData>();
	const string catalog_path = ResolveAnalysisCatalog(context, input, {"pvar", "psam"}, "plink_freq");
	if (catalog_path.empty()) {
		bind_data->pgen_path = input.inputs[0].GetValue<string>();
	}

	auto &fs = FileSystem::GetFileSystem(context);

	// --- Named parameters ---
	string build_str = "GRCh38"; // default genome build for PAR boundary detection
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "counts") {
			bind_data->include_counts = kv.second.GetValue<bool>();
		} else if (kv.first == "dosage") {
			bind_data->include_dosage = kv.second.GetValue<bool>();
		} else if (kv.first == "build") {
			build_str = kv.second.GetValue<string>();
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled after pgenlib init
		}
	}
	bind_data->par_bounds = ResolveParBounds(build_str, "plink_freq");

	if (!catalog_path.empty()) {
		BindFreqCatalog(context, input, *bind_data, catalog_path);
		SetFreqColumns(*bind_data, return_types, names);
		return std::move(bind_data);
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		bind_data->pvar_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_freq: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".psam", ".fam"});
		// .psam is optional for plink_freq — frequency computation only needs sample count
	}

	// --- Initialize pgenlib (Phase 1) to get counts ---
	LocalizePgenIfRequested(context, bind_data->pgen_path, bind_data->localize_guard);
	bind_data->use_vfs = PgenIoUseVfs(context, bind_data->pgen_path);
	PgenVfsScope pgen_vfs_scope(context, bind_data->use_vfs, PgenVfsRole::HEADER); // bind: the scan reuses its blocks

	// The decoded header (cache tier 1): shared with the scan threads and later queries.
	auto bind_header = AcquirePgenHeader(context, bind_data->use_vfs, bind_data->pgen_path, UINT32_MAX, UINT32_MAX,
	                                     "plink_freq: failed to open");
	const plink2::PgenFileInfo &pgfi = bind_header->pgfi;
	bind_data->raw_variant_ct = pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = pgfi.raw_sample_ct;
	bind_data->file_has_dosage = (pgfi.gflags & plink2::kfPgenGlobalDosagePresent) != 0;

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "plink_freq");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_freq: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// ParseSelection captures samples/region symbolically; ResolveSelection (below) composes
	// ResolveSampleIndices/BuildSampleSubset/ParseRegion over the loaded metadata.
	Selection selection = ParseSelection(input.named_parameters);

	// --- Load sample info (optional) ---
	// samples := [names] over a .psam.parquet reads only the named rows (IID, SEX), not every IID:
	// plink_freq uses no other sample's IID.
	if (!bind_data->psam_path.empty()) {
		vector<std::pair<string, uint32_t>> named_matches;
		const bool named_lookup =
		    selection.has_samples && LookupNamedSamplesInParquet(context, bind_data->psam_path, selection.samples,
		                                                         bind_data->sample_info, named_matches);
		if (!named_lookup) {
			bind_data->sample_info = LoadSampleMetadata(context, bind_data->psam_path);
		}
		bind_data->has_sample_info = true;

		if (static_cast<uint32_t>(bind_data->sample_info.sample_ct) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_freq: sample count mismatch: .pgen has %u samples, "
			                            ".psam/.fam '%s' has %llu samples",
			                            bind_data->raw_sample_ct, bind_data->psam_path,
			                            static_cast<unsigned long long>(bind_data->sample_info.sample_ct));
		}
		if (named_lookup) {
			selection.samples = NamedSampleIndices(selection.samples, named_matches, "plink_freq");
		}
	}

	// --- Resolve the axis selection (samples + region) via the shared handler ---
	// Samples resolve before region, psam is optional (sample_info may be null), and
	// count-mismatch checks remain plink_freq's own (above). BuildAlignedSex still runs here
	// with the resolved subset — its no-subset case must pass nullptr, not &subset_sorted.
	ResolvedSelection resolved =
	    ResolveSelection(selection, bind_data->variants, bind_data->has_sample_info ? &bind_data->sample_info : nullptr,
	                     bind_data->raw_sample_ct, bind_data->raw_variant_ct, "plink_freq");

	bind_data->effective_sample_ct = resolved.effective_sample_ct;
	if (resolved.has_sample_subset) {
		bind_data->sample_subset = std::move(resolved.sample_subset);
		bind_data->has_sample_subset = true;
	}
	bind_data->variant_range = resolved.variant_range;

	// --- Build ploidy-aware sex array (aligned to effective sample order) ---
	if (bind_data->has_sample_info && !bind_data->sample_info.sexes.empty()) {
		bind_data->aligned_sex =
		    BuildAlignedSex(bind_data->sample_info, resolved.has_subset ? &resolved.subset_sorted : nullptr);
		bind_data->have_sex = !bind_data->aligned_sex.empty();
	}

	SetFreqColumns(*bind_data, return_types, names);
	return std::move(bind_data);
}

static void SetFreqColumns(PlinkFreqBindData &bind_data_ref, vector<LogicalType> &return_types,
                           vector<CompatName> &names) {
	auto *bind_data = &bind_data_ref;
	names = {"CHROM", "POS", "ID", "REF", "ALT", "ALT_FREQ", "OBS_CT"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::INTEGER};

	if (bind_data->include_counts) {
		names.insert(names.end(), {"HOM_REF_CT", "HET_CT", "HOM_ALT_CT", "MISSING_CT"});
		return_types.insert(return_types.end(),
		                    {LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER});
	}

	if (bind_data->include_dosage) {
		bind_data->imp_r2_col_idx = names.size();
		names.push_back("IMP_R2");
		return_types.push_back(LogicalType::DOUBLE);
	}
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkFreqInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkFreqBindData>();
	auto state = make_uniq<PlinkFreqGlobalState>();

	if (bind_data.cohort) {
		// A catalog read scans positions over the cohort's segments.
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.cohort->scan_ct;
	} else if (bind_data.variant_range.has_filter) {
		state->start_variant_idx = bind_data.variant_range.start_idx;
		state->end_variant_idx = bind_data.variant_range.end_idx;
	} else {
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.raw_variant_ct;
	}

	state->next_variant_idx.store(state->start_variant_idx);
	state->column_ids = input.column_ids;
	state->max_threads_config = GetPlinkingMaxThreads(context);

	// Warn once per query, here rather than in the scan, because the scan runs on
	// several threads and only this point knows the final selected range.
	if (bind_data.cohort) {
		WarnOnMultiallelicCount(bind_data.cohort->CountMultiallelic(), "plink_freq");
	} else {
		WarnOnMultiallelicCollapse(bind_data.variants, state->start_variant_idx, state->end_variant_idx, "plink_freq");
	}

	// Check if any frequency/count/dosage columns are projected
	state->need_frequencies = false;
	for (auto col_id : input.column_ids) {
		if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}
		if (col_id >= COL_ALT_FREQ && col_id <= COL_MISSING_CT) {
			state->need_frequencies = true;
			break;
		}
		if (bind_data.include_dosage && col_id == bind_data.imp_r2_col_idx) {
			state->need_frequencies = true;
			break;
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkFreqInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                              GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkFreqBindData>();
	auto &gstate = global_state->Cast<PlinkFreqGlobalState>();
	auto state = make_uniq<PlinkFreqLocalState>();

	if (!gstate.need_frequencies) {
		return std::move(state);
	}

	if (bind_data.cohort) {
		// The scan opens one ConcatPgenReader per variant block (CatalogCohortReader). Only the
		// decode buffers for the sex-aware path are needed here, sized to the selected samples.
		uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(bind_data.effective_sample_ct);
		state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));
		std::memset(state->genovec_buf.ptr, 0, genovec_word_ct * sizeof(uintptr_t));
		state->geno_bytes.resize(bind_data.effective_sample_ct);
		state->cohort_ready = true;
		return std::move(state);
	}

	// --- Initialize per-thread PgenFileInfo + PgenReader ---
	PgenVfsScope pgen_vfs_scope(context.client, bind_data.use_vfs);
	plink2::PreinitPgr(&state->pgr);
	state->header = AcquirePgenHeader(context.client, bind_data.use_vfs, bind_data.pgen_path, bind_data.raw_variant_ct,
	                                  bind_data.raw_sample_ct, "plink_freq: thread init failed");
	InitPgenReaderOn(*state->header, bind_data.pgen_path, state->pgr, state->pgr_alloc_buf, "plink_freq");

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	// Buffers for decoding raw hardcalls on sex/organelle chromosomes.
	uintptr_t genovec_word_ct = plink2::NypCtToAlignedWordCt(bind_data.raw_sample_ct);
	state->genovec_buf.Allocate(genovec_word_ct * sizeof(uintptr_t));
	std::memset(state->genovec_buf.ptr, 0, genovec_word_ct * sizeof(uintptr_t));
	state->geno_bytes.resize(bind_data.effective_sample_ct);

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

		for (uint32_t pos = batch_start; pos < batch_end; pos++) {
			// One fileset: the scan position is the file's variant index. A catalog: the position
			// maps to a variant block and a block-local index, read through that block's reader.
			const VariantMetadataIndex *variants_ptr = &bind_data.variants;
			uint32_t vidx = pos;
			ConcatPgenReader *cat = nullptr;
			if (bind_data.cohort) {
				auto &seg = bind_data.cohort->SegmentAt(pos);
				variants_ptr = &bind_data.cohort->blocks[seg.block].variants;
				vidx = seg.start + (pos - seg.scan_start);
				if (lstate.cohort_ready) {
					cat = &lstate.cohort_reader.Open(context, *bind_data.cohort, seg.block, bind_data.cohort_has_subset,
					                                 bind_data.cohort_subset_sorted);
				}
			}
			const VariantMetadataIndex &variants = *variants_ptr;

			// Classify ploidy from the chromosome (and X position vs PAR). Autosomes
			// and the PAR keep the fast diploid PgrGetCounts path (no behavior change);
			// chrX non-PAR / chrY / chrMT take the ploidy/sex-aware hardcall path.
			ChromPloidy ploidy = ChromPloidy::AUTOSOMAL;
			if (gstate.need_frequencies) {
				ploidy = ClassifyChromPloidy(variants.GetChrom(vidx), variants.GetPos(vidx), bind_data.par_bounds);
			}
			bool sex_aware = (ploidy != ChromPloidy::AUTOSOMAL);

			// Compute genotype counts
			STD_ARRAY_DECL(uint32_t, 4, genocounts);
			genocounts[0] = genocounts[1] = genocounts[2] = genocounts[3] = 0;
			uint64_t all_dosages[2] = {0, 0};
			double imp_r2 = 0.0;
			SexAwareCounts sac;

			if (gstate.need_frequencies && (lstate.initialized || cat)) {
				if (sex_aware) {
					plink2::PglErr err = cat ? cat->Get(vidx, lstate.genovec_buf.As<uintptr_t>())
					                         : plink2::PgrGet(sample_include, lstate.pssi, sample_ct, vidx, &lstate.pgr,
					                                          lstate.genovec_buf.As<uintptr_t>());
					if (err != plink2::kPglRetSuccess) {
						throw IOException("plink_freq: PgrGet failed for variant %u", vidx);
					}
					plink2::GenoarrToBytesMinus9(lstate.genovec_buf.As<uintptr_t>(), sample_ct,
					                             lstate.geno_bytes.data());
					const uint8_t *sex_ptr = bind_data.have_sex ? bind_data.aligned_sex.data() : nullptr;
					sac =
					    ComputeSexAwareCounts(lstate.geno_bytes.data(), sample_ct, ploidy, sex_ptr, bind_data.have_sex);
				} else if (bind_data.include_dosage) {
					plink2::PglErr err =
					    plink2::PgrGetDCounts(sample_include, interleaved_vec, lstate.pssi, sample_ct, vidx,
					                          0, // is_minimac3_r2 = 0 (standard R²)
					                          &lstate.pgr, &imp_r2, genocounts, all_dosages);
					if (err != plink2::kPglRetSuccess) {
						throw IOException("plink_freq: PgrGetDCounts failed for variant %u", vidx);
					}
				} else {
					plink2::PglErr err = cat ? cat->GetCounts(vidx, genocounts)
					                         : plink2::PgrGetCounts(sample_include, interleaved_vec, lstate.pssi,
					                                                sample_ct, vidx, &lstate.pgr, genocounts);
					if (err != plink2::kPglRetSuccess) {
						throw IOException("plink_freq: PgrGetCounts failed for variant %u", vidx);
					}
				}
			}

			// Compute derived values
			// kDosageMid = 16384: each non-missing sample contributes 2*kDosageMid to
			// all_dosages total (REF + ALT always sum to 2.0 per sample)
			static constexpr uint64_t kDosageMid = 16384;

			uint32_t hardcall_obs_sample_ct = genocounts[0] + genocounts[1] + genocounts[2];
			uint32_t obs_ct;
			double alt_freq;
			bool freq_is_null = false;

			// Count-column sources (may differ from genocounts on sex chromosomes).
			int32_t out_hom_ref = static_cast<int32_t>(genocounts[0]);
			int32_t out_het = static_cast<int32_t>(genocounts[1]);
			int32_t out_hom_alt = static_cast<int32_t>(genocounts[2]);
			int32_t out_missing = static_cast<int32_t>(genocounts[3]);
			bool counts_are_null = false;

			if (sex_aware) {
				// Ploidy/sex-aware allele frequency: haploid samples contribute 1
				// allele; OBS_CT is an allele count (not 2×sample_ct).
				if (sac.sex_unavailable || sac.obs_allele_ct == 0) {
					freq_is_null = true;
					alt_freq = 0.0;
					obs_ct = 0;
					counts_are_null = sac.sex_unavailable; // no sex → genotype counts undefined too
				} else {
					alt_freq = static_cast<double>(sac.alt_allele_ct) / static_cast<double>(sac.obs_allele_ct);
					obs_ct = sac.obs_allele_ct;
				}
				out_hom_ref = static_cast<int32_t>(sac.geno_hom_ref);
				out_het = static_cast<int32_t>(sac.geno_het);
				out_hom_alt = static_cast<int32_t>(sac.geno_hom_alt);
				out_missing = static_cast<int32_t>(sac.geno_missing);
			} else if (bind_data.include_dosage) {
				// Dosage-weighted: all_dosages values are scaled by kDosageMid
				uint64_t total_dosage = all_dosages[0] + all_dosages[1];
				if (total_dosage == 0) {
					freq_is_null = true;
					alt_freq = 0.0;
					obs_ct = 0;
				} else {
					alt_freq = static_cast<double>(all_dosages[1]) / static_cast<double>(total_dosage);
					// Effective allele observations: total_dosage / kDosageMid
					// Always an integer since each sample contributes exactly 2*kDosageMid
					obs_ct = static_cast<uint32_t>(total_dosage / kDosageMid);
				}
			} else if (hardcall_obs_sample_ct == 0) {
				freq_is_null = true;
				alt_freq = 0.0;
				obs_ct = 0;
			} else {
				obs_ct = 2 * hardcall_obs_sample_ct;
				alt_freq = (static_cast<double>(genocounts[1]) + 2.0 * static_cast<double>(genocounts[2])) /
				           (2.0 * static_cast<double>(hardcall_obs_sample_ct));
			}

			// Fill projected columns
			for (idx_t out_col = 0; out_col < column_ids.size(); out_col++) {
				auto file_col = column_ids[out_col];
				if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
					continue;
				}

				auto &vec = output.data[out_col];

				// Handle dynamic IMP_R2 column before the switch to avoid collision
				// with static column constants (imp_r2_col_idx may equal COL_HOM_REF_CT
				// when counts is disabled)
				if (bind_data.include_dosage && file_col == bind_data.imp_r2_col_idx) {
					if (sex_aware) {
						// IMP_R2 is not computed on the ploidy-aware hardcall path.
						FlatVector::SetNull(vec, rows_emitted, true);
					} else if (bind_data.file_has_dosage) {
						// Note: file-level check. In files with mixed dosage/hardcall variants,
						// PgrGetDCounts returns a hardcall-derived R² approximation for variants
						// without dosage entries. Real imputed files typically have dosage for all
						// variants, so this covers the common cases.
						CompatFlatDataMutable<double>(vec)[rows_emitted] = imp_r2;
					} else {
						// No dosage data in file — IMP_R2 is not meaningful
						FlatVector::SetNull(vec, rows_emitted, true);
					}
					continue;
				}

				switch (file_col) {
				case COL_CHROM: {
					auto val = variants.GetChrom(vidx);
					CompatFlatDataMutable<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case COL_POS: {
					CompatFlatDataMutable<int32_t>(vec)[rows_emitted] = variants.GetPos(vidx);
					break;
				}
				case COL_ID: {
					auto val = variants.GetId(vidx);
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_REF: {
					auto val = variants.GetRef(vidx);
					CompatFlatDataMutable<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case COL_ALT: {
					auto val = variants.GetAlt(vidx);
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_ALT_FREQ: {
					if (freq_is_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<double>(vec)[rows_emitted] = alt_freq;
					}
					break;
				}
				case COL_OBS_CT: {
					CompatFlatDataMutable<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(obs_ct);
					break;
				}
				case COL_HOM_REF_CT: {
					if (counts_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<int32_t>(vec)[rows_emitted] = out_hom_ref;
					}
					break;
				}
				case COL_HET_CT: {
					if (counts_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<int32_t>(vec)[rows_emitted] = out_het;
					}
					break;
				}
				case COL_HOM_ALT_CT: {
					if (counts_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<int32_t>(vec)[rows_emitted] = out_hom_alt;
					}
					break;
				}
				case COL_MISSING_CT: {
					if (counts_are_null) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						CompatFlatDataMutable<int32_t>(vec)[rows_emitted] = out_missing;
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

	CompatSetOutputCardinality(output, rows_emitted);
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

void RegisterPlinkFreq(ExtensionLoader &loader) {
	// Two overloads: a .pgen or catalog path, and none (catalog := '<path>' names the inputs).
	TableFunctionSet set("plink_freq");
	for (auto &args : vector<vector<LogicalType>> {{LogicalType::VARCHAR}, {}}) {
		TableFunction plink_freq("plink_freq", args, PlinkFreqScan, PlinkFreqBind, PlinkFreqInitGlobal,
		                         PlinkFreqInitLocal);

		plink_freq.projection_pushdown = true;

		DeclareNamedParameters(plink_freq, {{"pvar", LogicalType::VARCHAR},
		                                    {"psam", LogicalType::VARCHAR},
		                                    {"samples", LogicalType::ANY},
		                                    {"region", LogicalType::VARCHAR},
		                                    {"selection", LogicalType::ANY},
		                                    {"counts", LogicalType::BOOLEAN},
		                                    {"dosage", LogicalType::BOOLEAN},
		                                    {"build", LogicalType::VARCHAR},
		                                    {"catalog", LogicalType::ANY}});
		CompatSetCaptureArgumentAliases(plink_freq);
		set.AddFunction(plink_freq);
	}

	RegisterTableSetWithDesc(loader, set, {}, "Compute allele frequencies and counts from PLINK files.",
	                         {"SELECT * FROM plink_freq('data/cohort.pgen')"});
}

} // namespace duckdb
