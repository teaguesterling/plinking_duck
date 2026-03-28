#include "plink_glm.hpp"
#include "plink_common.hpp"
#include "psam_reader.hpp"
#include "plink2_glm_logistic_math.hpp"

#include <atomic>
#include <cmath>
#include <cstring>

namespace duckdb {

// ---------------------------------------------------------------------------
// Psam column extraction helper
// ---------------------------------------------------------------------------

//! Extract a named column from pre-read .psam/.fam lines as a vector of doubles.
//! Returns one value per sample in file order. "NA" and missing values
//! become NAN. Throws on column not found or non-numeric values.
static vector<double> LoadPsamColumnAsDouble(const vector<string> &lines, const PsamHeaderInfo &header,
                                             const string &psam_path, const string &column_name) {

	// Find the column index
	idx_t col_idx = DConstants::INVALID_INDEX;
	for (idx_t i = 0; i < header.column_names.size(); i++) {
		if (header.column_names[i] == column_name) {
			col_idx = i;
			break;
		}
	}
	if (col_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("plink_glm: psam file '%s' has no column '%s'", psam_path, column_name);
	}

	// Data lines start at index 1 for .psam (skip header), 0 for .fam
	idx_t data_start = (header.format != PsamFormat::FAM) ? 1 : 0;

	vector<double> values;
	for (idx_t i = data_start; i < lines.size(); i++) {
		auto &line = lines[i];
		if (line.empty()) {
			continue;
		}

		// Split by tab (psam) or whitespace (fam)
		vector<string> fields;
		if (header.format == PsamFormat::FAM) {
			// Space-delimited
			idx_t pos = 0;
			while (pos < line.size()) {
				while (pos < line.size() && (line[pos] == ' ' || line[pos] == '\t')) {
					pos++;
				}
				if (pos >= line.size()) {
					break;
				}
				idx_t start = pos;
				while (pos < line.size() && line[pos] != ' ' && line[pos] != '\t') {
					pos++;
				}
				fields.push_back(line.substr(start, pos - start));
			}
		} else {
			// Tab-delimited
			idx_t start = 0;
			for (idx_t j = 0; j <= line.size(); j++) {
				if (j == line.size() || line[j] == '\t') {
					fields.push_back(line.substr(start, j - start));
					start = j + 1;
				}
			}
		}

		if (col_idx >= fields.size()) {
			throw IOException("plink_glm: psam file '%s' line %d has %d fields, expected at least %d", psam_path, i + 1,
			                  fields.size(), col_idx + 1);
		}

		auto &field = fields[col_idx];
		if (field == "NA" || field == "na" || field == "-9") {
			values.push_back(NAN);
		} else {
			try {
				values.push_back(std::stod(field));
			} catch (const std::exception &) {
				throw InvalidInputException("plink_glm: psam column '%s' has non-numeric value '%s' at line %d",
				                            column_name, field, i + 1);
			}
		}
	}

	return values;
}

// ---------------------------------------------------------------------------
// Column indices
// ---------------------------------------------------------------------------

static constexpr idx_t COL_CHROM = 0;
static constexpr idx_t COL_POS = 1;
static constexpr idx_t COL_ID = 2;
static constexpr idx_t COL_REF = 3;
static constexpr idx_t COL_ALT = 4;
static constexpr idx_t COL_A1 = 5;
static constexpr idx_t COL_A1_FREQ = 6;
static constexpr idx_t COL_TEST = 7;
static constexpr idx_t COL_OBS_CT = 8;
static constexpr idx_t COL_BETA = 9;
static constexpr idx_t COL_SE = 10;
static constexpr idx_t COL_T_STAT = 11;
static constexpr idx_t COL_P = 12;
static constexpr idx_t COL_ERRCODE = 13;
static constexpr idx_t COL_OR = 14;
static constexpr idx_t COL_FIRTH_YN = 15;

// ---------------------------------------------------------------------------
// T-distribution p-value via regularized incomplete beta function
// ---------------------------------------------------------------------------

static double LogGamma(double x) {
	static const double c[] = {0.99999999999980993,  676.5203681218851,     -1259.1392167224028,
	                           771.32342877765313,   -176.61502916214059,   12.507343278686905,
	                           -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
	if (x < 0.5) {
		return std::log(M_PI / std::sin(M_PI * x)) - LogGamma(1.0 - x);
	}
	x -= 1.0;
	double a = c[0];
	double t = x + 7.5;
	for (int i = 1; i < 9; i++) {
		a += c[i] / (x + i);
	}
	return 0.5 * std::log(2.0 * M_PI) + (x + 0.5) * std::log(t) - t + std::log(a);
}

//! Regularized incomplete beta function I_x(a, b) via continued fraction.
static double BetaIncomplete(double a, double b, double x) {
	if (x <= 0.0) {
		return 0.0;
	}
	if (x >= 1.0) {
		return 1.0;
	}

	// Use symmetry relation for better convergence
	if (x > (a + 1.0) / (a + b + 2.0)) {
		return 1.0 - BetaIncomplete(b, a, 1.0 - x);
	}

	double lbeta = LogGamma(a) + LogGamma(b) - LogGamma(a + b);
	double prefix = std::exp(a * std::log(x) + b * std::log(1.0 - x) - lbeta) / a;

	static constexpr double TINY = 1e-30;
	static constexpr double EPS = 1e-14;
	static constexpr int MAX_ITER = 200;

	double c_val = 1.0;
	double d = 1.0 - (a + b) * x / (a + 1.0);
	if (std::abs(d) < TINY) {
		d = TINY;
	}
	d = 1.0 / d;
	double result = d;

	for (int m = 1; m <= MAX_ITER; m++) {
		// Even step
		double num = m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
		d = 1.0 + num * d;
		if (std::abs(d) < TINY) {
			d = TINY;
		}
		c_val = 1.0 + num / c_val;
		if (std::abs(c_val) < TINY) {
			c_val = TINY;
		}
		d = 1.0 / d;
		result *= d * c_val;

		// Odd step
		num = -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
		d = 1.0 + num * d;
		if (std::abs(d) < TINY) {
			d = TINY;
		}
		c_val = 1.0 + num / c_val;
		if (std::abs(c_val) < TINY) {
			c_val = TINY;
		}
		d = 1.0 / d;
		double delta = d * c_val;
		result *= delta;

		if (std::abs(delta - 1.0) < EPS) {
			break;
		}
	}

	return prefix * result;
}

//! Two-tailed p-value from t-statistic with df degrees of freedom.
//! Uses the identity: P = I_x(df/2, 1/2) where x = df/(df + t^2).
static double TstatToPvalue(double t_stat, double df) {
	if (df <= 0.0 || !std::isfinite(t_stat)) {
		return NAN;
	}
	if (t_stat == 0.0) {
		return 1.0;
	}
	double x = df / (df + t_stat * t_stat);
	return BetaIncomplete(df / 2.0, 0.5, x);
}

// ---------------------------------------------------------------------------
// Normal distribution p-value (for logistic regression z-statistic)
// ---------------------------------------------------------------------------

//! Standard normal CDF via error function.
static double NormalCdf(double x) {
	return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

//! Two-tailed p-value from z-statistic (standard normal).
static double ZstatToPvalue(double z_stat) {
	if (!std::isfinite(z_stat)) {
		return NAN;
	}
	if (z_stat == 0.0) {
		return 1.0;
	}
	return 2.0 * (1.0 - NormalCdf(std::abs(z_stat)));
}

// plink2's LinearRegressionInv is a header-inline that calls into plink2_matrix.
// Under NOLAPACK, MatrixInvertBuf1 is a typedef in global scope (not in plink2::).
using plink2::InvertSymmdefFmatrixFirstHalf;
using plink2::InvertSymmdefFmatrixSecondHalf;
using plink2::kFloatPerFVec;
using plink2::LinearRegressionInv;
using plink2::MultiplySelfTranspose;
using plink2::ReflectFmatrix0;
using plink2::RoundUpPow2;

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkGlmBindData : public TableFunctionData {
	string pgen_path;
	string pvar_path;
	string psam_path;

	VariantMetadataIndex variants;

	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Sample subsetting
	bool has_sample_subset = false;
	unique_ptr<SampleSubset> sample_subset;
	uint32_t effective_sample_ct = 0;

	// Region filtering
	VariantRange variant_range;

	// Phenotype values aligned with effective samples (NaN = missing)
	vector<double> phenotype;

	// Covariates: covariate_values[covar_idx][sample_idx] for effective samples
	vector<string> covariate_names;
	vector<vector<double>> covariate_values;

	// Model
	bool is_logistic = false;
	bool use_firth = true;
	uint32_t predictor_ct = 0; // intercept + genotype + covariates

	// P-value threshold filter (NaN = no filter)
	double p_threshold = std::numeric_limits<double>::quiet_NaN();
};

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------

struct PlinkGlmGlobalState : public GlobalTableFunctionState {
	std::atomic<uint32_t> next_variant_idx {0};
	uint32_t start_variant_idx = 0;
	uint32_t end_variant_idx = 0;
	vector<column_t> column_ids;
	bool need_regression = false;
	uint32_t max_threads_config = 0;

	idx_t MaxThreads() const override {
		uint32_t range = end_variant_idx - start_variant_idx;
		return ApplyMaxThreadsCap(range / 500 + 1, max_threads_config);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

// Pre-allocated scratch buffers for logistic regression, reused across variants.
// Allocated once per thread in InitLocal based on max sample count and predictor count.
struct LogisticBuffers {
	// Core IRLS buffers
	vector<float> xx;    // p * max_sample_ctav (covariate-major design matrix)
	vector<float> yy;    // max_sample_ctav (phenotype)
	vector<float> coef;  // predictor_ctav
	vector<float> ll;    // p * predictor_ctav (Cholesky output)
	vector<float> pp;    // max_sample_ctav (predicted probabilities)
	vector<float> vv;    // max_sample_ctav (variances)
	vector<float> hh;    // p * predictor_ctav (Hessian / inverted var-cov)
	vector<float> grad;  // predictor_ctav
	vector<float> dcoef; // predictor_ctav

	// Firth-specific buffers
	vector<double> half_inverted_buf;    // p * p
	vector<MatrixInvertBuf1> inv_1d_buf; // 2 * p
	vector<double> dbl_2d_buf;           // p * p
	vector<float> ustar;                 // predictor_ctav
	vector<float> delta;                 // predictor_ctav
	vector<float> hdiag;                 // max_sample_ctav
	vector<float> ww;                    // max_sample_ctav
	vector<float> hh0;                   // p * predictor_ctav
	vector<float> tmpnxk_buf;            // p * max_sample_ctav

	// Hessian inversion buffers (reused for non-Firth and Firth paths)
	// half_inverted_buf, inv_1d_buf, dbl_2d_buf are shared above

	// Non-missing sample indices scratch
	vector<uint32_t> nm_indices;

	bool allocated = false;

	void Allocate(uint32_t max_sample_ct, uint32_t predictor_ct, bool use_firth) {
		uint32_t max_sample_ctav = RoundUpPow2(max_sample_ct, kFloatPerFVec);
		uint32_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
		uint32_t p = predictor_ct;

		xx.resize(p * max_sample_ctav);
		yy.resize(max_sample_ctav);
		coef.resize(predictor_ctav);
		ll.resize(p * predictor_ctav);
		pp.resize(max_sample_ctav);
		vv.resize(max_sample_ctav);
		hh.resize(p * predictor_ctav);
		grad.resize(predictor_ctav);
		dcoef.resize(predictor_ctav);

		// Hessian inversion buffers (needed for all logistic paths)
		half_inverted_buf.resize(p * p);
		inv_1d_buf.resize(2 * p);
		dbl_2d_buf.resize(p * p);

		if (use_firth) {
			ustar.resize(predictor_ctav);
			delta.resize(predictor_ctav);
			hdiag.resize(max_sample_ctav);
			ww.resize(max_sample_ctav);
			hh0.resize(p * predictor_ctav);
			tmpnxk_buf.resize(p * max_sample_ctav);
		}

		nm_indices.resize(max_sample_ct);

		allocated = true;
	}

	//! Zero all buffers before reuse. Called once per variant.
	void Reset(uint32_t n, uint32_t p, bool use_firth) {
		uint32_t sample_ctav = RoundUpPow2(n, kFloatPerFVec);
		uint32_t predictor_ctav = RoundUpPow2(p, kFloatPerFVec);

		std::fill_n(xx.data(), p * sample_ctav, 0.0f);
		std::fill_n(yy.data(), sample_ctav, 0.0f);
		std::fill_n(coef.data(), predictor_ctav, 0.0f);
		std::fill_n(ll.data(), p * predictor_ctav, 0.0f);
		std::fill_n(pp.data(), sample_ctav, 0.0f);
		std::fill_n(vv.data(), sample_ctav, 0.0f);
		std::fill_n(hh.data(), p * predictor_ctav, 0.0f);
		std::fill_n(grad.data(), predictor_ctav, 0.0f);
		std::fill_n(dcoef.data(), predictor_ctav, 0.0f);
	}
};

struct PlinkGlmLocalState : public LocalTableFunctionState {
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

	// Pre-allocated logistic regression scratch buffers
	LogisticBuffers logistic_bufs;

	bool initialized = false;

	~PlinkGlmLocalState() {
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

static unique_ptr<FunctionData> PlinkGlmBind(ClientContext &context, TableFunctionBindInput &input,
                                             vector<LogicalType> &return_types, vector<string> &names) {
	auto bind_data = make_uniq<PlinkGlmBindData>();
	string prefix = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);

	string model_str = "auto";

	// --- Named parameters (first pass: file paths and options) ---
	for (auto &kv : input.named_parameters) {
		if (kv.first == "pvar") {
			bind_data->pvar_path = kv.second.GetValue<string>();
		} else if (kv.first == "psam") {
			bind_data->psam_path = kv.second.GetValue<string>();
		} else if (kv.first == "model") {
			model_str = kv.second.GetValue<string>();
		} else if (kv.first == "firth") {
			bind_data->use_firth = kv.second.GetValue<bool>();
		} else if (kv.first == "p_threshold") {
			double pt = kv.second.GetValue<double>();
			if (pt <= 0.0 || pt > 1.0) {
				throw InvalidInputException("plink_glm: p_threshold must be in (0, 1], got %g", pt);
			}
			bind_data->p_threshold = pt;
		}
	}

	// --- Resolve pgen path from prefix ---
	{
		string candidate = prefix + ".pgen";
		if (fs.FileExists(candidate)) {
			bind_data->pgen_path = candidate;
		} else if (fs.FileExists(prefix)) {
			bind_data->pgen_path = prefix;
		} else {
			throw InvalidInputException("plink_glm: cannot find .pgen file for prefix '%s' (tried '%s')", prefix,
			                            candidate);
		}
	}

	// --- Auto-discover companion files ---
	if (bind_data->pvar_path.empty()) {
		// Try prefix-based discovery first, then companion file discovery
		for (auto &ext : {".pvar", ".bim"}) {
			auto candidate = prefix + ext;
			if (fs.FileExists(candidate)) {
				bind_data->pvar_path = candidate;
				break;
			}
		}
		if (bind_data->pvar_path.empty()) {
			bind_data->pvar_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".pvar", ".bim"});
		}
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_glm: cannot find .pvar or .bim file for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            prefix);
		}
	}

	if (bind_data->psam_path.empty()) {
		for (auto &ext : {".psam", ".fam"}) {
			auto candidate = prefix + ext;
			if (fs.FileExists(candidate)) {
				bind_data->psam_path = candidate;
				break;
			}
		}
		if (bind_data->psam_path.empty()) {
			bind_data->psam_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".psam", ".fam"});
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
		throw IOException("plink_glm: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
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
		throw IOException("plink_glm: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}

	// --- Load variant metadata ---
	bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "plink_glm");

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_glm: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Process samples parameter ---
	bind_data->effective_sample_ct = bind_data->raw_sample_ct;

	auto samples_it = input.named_parameters.find("samples");
	if (samples_it != input.named_parameters.end()) {
		if (bind_data->psam_path.empty()) {
			throw InvalidInputException("plink_glm: samples parameter requires .psam/.fam file");
		}
		auto sample_info = LoadSampleMetadata(context, bind_data->psam_path);
		auto indices = ResolveSampleIndices(samples_it->second, bind_data->raw_sample_ct, &sample_info, "plink_glm");
		bind_data->sample_subset = make_uniq<SampleSubset>(BuildSampleSubset(bind_data->raw_sample_ct, indices));
		bind_data->has_sample_subset = true;
		bind_data->effective_sample_ct = bind_data->sample_subset->subset_sample_ct;
	}

	// --- Process region parameter ---
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		bind_data->variant_range = ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_glm");
	}

	// --- Process phenotype parameter ---
	auto pheno_it = input.named_parameters.find("phenotype");
	if (pheno_it == input.named_parameters.end()) {
		throw InvalidInputException("plink_glm: phenotype parameter is required");
	}

	auto &pheno_val = pheno_it->second;

	// Pre-read psam file once if any parameter uses column name dispatch
	vector<string> psam_lines;
	PsamHeaderInfo psam_header;
	bool psam_loaded = false;

	auto EnsurePsamLoaded = [&]() {
		if (!psam_loaded && !bind_data->psam_path.empty()) {
			psam_lines = ReadFileLines(context, bind_data->psam_path);
			psam_header = ParsePsamHeaderFromLines(psam_lines, bind_data->psam_path);
			psam_loaded = true;
		}
	};

	// Type dispatch: VARCHAR → psam column name, LIST → direct values
	if (pheno_val.type().id() == LogicalTypeId::VARCHAR) {
		// --- Phenotype from psam column ---
		if (bind_data->psam_path.empty()) {
			throw InvalidInputException("plink_glm: phenotype as column name requires a .psam/.fam file");
		}
		EnsurePsamLoaded();
		string pheno_col = pheno_val.GetValue<string>();
		auto raw_values = LoadPsamColumnAsDouble(psam_lines, psam_header, bind_data->psam_path, pheno_col);
		if (static_cast<uint32_t>(raw_values.size()) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_glm: psam has %llu samples but .pgen has %u",
			                            static_cast<unsigned long long>(raw_values.size()), bind_data->raw_sample_ct);
		}

		bind_data->phenotype.resize(bind_data->effective_sample_ct);
		if (bind_data->has_sample_subset) {
			const uintptr_t *sample_include = bind_data->sample_subset->SampleInclude();
			uint32_t out_idx = 0;
			for (uint32_t raw_idx = 0; raw_idx < bind_data->raw_sample_ct; raw_idx++) {
				if (plink2::IsSet(sample_include, raw_idx)) {
					bind_data->phenotype[out_idx] = raw_values[raw_idx];
					out_idx++;
				}
			}
		} else {
			bind_data->phenotype = std::move(raw_values);
		}
	} else if (pheno_val.type().id() == LogicalTypeId::LIST) {
		// --- Phenotype from direct LIST(DOUBLE) ---
		auto &pheno_children = ListValue::GetChildren(pheno_val);

		if (static_cast<uint32_t>(pheno_children.size()) != bind_data->raw_sample_ct) {
			throw InvalidInputException("plink_glm: phenotype length (%llu) must match sample count (%u)",
			                            static_cast<unsigned long long>(pheno_children.size()),
			                            bind_data->raw_sample_ct);
		}

		bind_data->phenotype.resize(bind_data->effective_sample_ct);
		if (bind_data->has_sample_subset) {
			const uintptr_t *sample_include = bind_data->sample_subset->SampleInclude();
			uint32_t out_idx = 0;
			for (uint32_t raw_idx = 0; raw_idx < bind_data->raw_sample_ct; raw_idx++) {
				if (plink2::IsSet(sample_include, raw_idx)) {
					if (pheno_children[raw_idx].IsNull()) {
						bind_data->phenotype[out_idx] = NAN;
					} else {
						bind_data->phenotype[out_idx] = pheno_children[raw_idx].GetValue<double>();
					}
					out_idx++;
				}
			}
		} else {
			for (uint32_t i = 0; i < bind_data->raw_sample_ct; i++) {
				if (pheno_children[i].IsNull()) {
					bind_data->phenotype[i] = NAN;
				} else {
					bind_data->phenotype[i] = pheno_children[i].GetValue<double>();
				}
			}
		}
	} else {
		throw InvalidInputException("plink_glm: phenotype must be a LIST(DOUBLE) or a VARCHAR column name");
	}

	// Validate phenotype
	uint32_t non_missing_ct = 0;
	double min_val = std::numeric_limits<double>::infinity();
	double max_val = -std::numeric_limits<double>::infinity();

	for (uint32_t i = 0; i < bind_data->effective_sample_ct; i++) {
		double v = bind_data->phenotype[i];
		if (!std::isnan(v)) {
			min_val = std::min(min_val, v);
			max_val = std::max(max_val, v);
			non_missing_ct++;
		}
	}

	if (non_missing_ct < 3) {
		throw InvalidInputException("plink_glm: need at least 3 non-missing phenotype values, got %u", non_missing_ct);
	}

	if (min_val == max_val) {
		throw InvalidInputException("plink_glm: constant phenotype (all values are %g)", min_val);
	}

	// --- Process covariates parameter ---
	auto covar_it = input.named_parameters.find("covariates");
	if (covar_it != input.named_parameters.end()) {
		auto &covar_val = covar_it->second;
		auto &covar_type = covar_val.type();

		if (covar_type.id() == LogicalTypeId::LIST &&
		    ListType::GetChildType(covar_type).id() == LogicalTypeId::VARCHAR) {
			// --- Covariates from psam column names: LIST(VARCHAR) ---
			if (bind_data->psam_path.empty()) {
				throw InvalidInputException("plink_glm: covariates as column names requires a .psam/.fam file");
			}
			EnsurePsamLoaded();

			auto &list_children = ListValue::GetChildren(covar_val);
			for (idx_t ci = 0; ci < list_children.size(); ci++) {
				string col_name = list_children[ci].GetValue<string>();
				auto raw_values = LoadPsamColumnAsDouble(psam_lines, psam_header, bind_data->psam_path, col_name);
				if (static_cast<uint32_t>(raw_values.size()) != bind_data->raw_sample_ct) {
					throw InvalidInputException("plink_glm: psam has %llu samples but .pgen has %u",
					                            static_cast<unsigned long long>(raw_values.size()),
					                            bind_data->raw_sample_ct);
				}

				// Check for NAN (which LoadPsamColumnAsDouble uses for NA/-9)
				for (uint32_t i = 0; i < raw_values.size(); i++) {
					if (std::isnan(raw_values[i])) {
						throw InvalidInputException("plink_glm: covariate '%s' contains NULL at sample %u", col_name,
						                            i);
					}
				}

				// Apply sample subset
				vector<double> covar_vals(bind_data->effective_sample_ct);
				if (bind_data->has_sample_subset) {
					const uintptr_t *si = bind_data->sample_subset->SampleInclude();
					uint32_t out_idx = 0;
					for (uint32_t raw_idx = 0; raw_idx < bind_data->raw_sample_ct; raw_idx++) {
						if (plink2::IsSet(si, raw_idx)) {
							covar_vals[out_idx] = raw_values[raw_idx];
							out_idx++;
						}
					}
				} else {
					covar_vals = std::move(raw_values);
				}

				bind_data->covariate_names.push_back(col_name);
				bind_data->covariate_values.push_back(std::move(covar_vals));
			}
		} else if (covar_type.id() == LogicalTypeId::STRUCT) {
			// --- Covariates from STRUCT of LIST(DOUBLE) ---
			auto &struct_children = StructType::GetChildTypes(covar_type);
			auto &struct_vals = StructValue::GetChildren(covar_val);

			for (idx_t ci = 0; ci < struct_children.size(); ci++) {
				auto &name = struct_children[ci].first;
				auto &val = struct_vals[ci];

				if (val.type().id() != LogicalTypeId::LIST) {
					throw InvalidInputException("plink_glm: covariate '%s' must be a LIST, got %s", name,
					                            val.type().ToString());
				}

				auto &list_children = ListValue::GetChildren(val);
				if (static_cast<uint32_t>(list_children.size()) != bind_data->raw_sample_ct) {
					throw InvalidInputException("plink_glm: covariate '%s' length (%llu) must match sample count (%u)",
					                            name, static_cast<unsigned long long>(list_children.size()),
					                            bind_data->raw_sample_ct);
				}

				// Extract values for effective samples (NULLs not allowed in covariates)
				vector<double> covar_vals(bind_data->effective_sample_ct);
				if (bind_data->has_sample_subset) {
					const uintptr_t *si = bind_data->sample_subset->SampleInclude();
					uint32_t out_idx = 0;
					for (uint32_t raw_idx = 0; raw_idx < bind_data->raw_sample_ct; raw_idx++) {
						if (plink2::IsSet(si, raw_idx)) {
							if (list_children[raw_idx].IsNull()) {
								throw InvalidInputException("plink_glm: covariate '%s' contains NULL at index %u", name,
								                            raw_idx);
							}
							covar_vals[out_idx] = list_children[raw_idx].GetValue<double>();
							out_idx++;
						}
					}
				} else {
					for (uint32_t i = 0; i < bind_data->raw_sample_ct; i++) {
						if (list_children[i].IsNull()) {
							throw InvalidInputException("plink_glm: covariate '%s' contains NULL at index %u", name, i);
						}
						covar_vals[i] = list_children[i].GetValue<double>();
					}
				}

				bind_data->covariate_names.push_back(name);
				bind_data->covariate_values.push_back(std::move(covar_vals));
			}
		} else {
			throw InvalidInputException("plink_glm: covariates must be a STRUCT of named LIST(DOUBLE) "
			                            "columns (e.g. {'age': list(age)}), or a LIST(VARCHAR) of "
			                            "psam column names (e.g. ['age', 'sex'])");
		}
	}

	// --- Model detection ---
	if (model_str == "linear") {
		bind_data->is_logistic = false;
	} else if (model_str == "logistic") {
		bind_data->is_logistic = true;
	} else if (model_str == "auto") {
		// Auto-detect: binary if all non-missing values are 0/1 or 1/2
		bool all_binary_01 = true;
		bool all_binary_12 = true;
		for (uint32_t i = 0; i < bind_data->effective_sample_ct; i++) {
			double v = bind_data->phenotype[i];
			if (std::isnan(v)) {
				continue;
			}
			if (v != 0.0 && v != 1.0) {
				all_binary_01 = false;
			}
			if (v != 1.0 && v != 2.0) {
				all_binary_12 = false;
			}
		}

		if (all_binary_01) {
			bind_data->is_logistic = true;
		} else if (all_binary_12) {
			// Recode 1/2 → 0/1
			bind_data->is_logistic = true;
			for (uint32_t i = 0; i < bind_data->effective_sample_ct; i++) {
				if (!std::isnan(bind_data->phenotype[i])) {
					bind_data->phenotype[i] -= 1.0;
				}
			}
		}
	} else {
		throw InvalidInputException("plink_glm: model must be 'auto', 'linear', or 'logistic', got '%s'", model_str);
	}

	bind_data->predictor_ct = static_cast<uint32_t>(2 + bind_data->covariate_values.size());

	// --- Output columns ---
	names = {"CHROM",  "POS",  "ID", "REF",    "ALT", "A1",      "A1_FREQ", "TEST",
	         "OBS_CT", "BETA", "SE", "T_STAT", "P",   "ERRCODE", "OR",      "FIRTH_YN"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::VARCHAR,
	                LogicalType::INTEGER, LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,
	                LogicalType::DOUBLE,  LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::VARCHAR};

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkGlmInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkGlmBindData>();
	auto state = make_uniq<PlinkGlmGlobalState>();

	if (bind_data.variant_range.has_filter) {
		state->start_variant_idx = bind_data.variant_range.start_idx;
		state->end_variant_idx = bind_data.variant_range.end_idx;
	} else {
		state->start_variant_idx = 0;
		state->end_variant_idx = bind_data.raw_variant_ct;
	}

	state->next_variant_idx.store(state->start_variant_idx);
	state->column_ids = input.column_ids;
	state->max_threads_config = GetPlinkingMaxThreads(context);

	// Check if any regression columns are projected, or if p_threshold requires it
	state->need_regression = !std::isnan(bind_data.p_threshold);
	for (auto col_id : input.column_ids) {
		if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}
		if (col_id >= COL_A1_FREQ) {
			state->need_regression = true;
			break;
		}
	}

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkGlmInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                             GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkGlmBindData>();
	auto &gstate = global_state->Cast<PlinkGlmGlobalState>();
	auto state = make_uniq<PlinkGlmLocalState>();

	if (!gstate.need_regression) {
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
		throw IOException("plink_glm: thread init failed (phase 1): %s", errstr_buf);
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
		throw IOException("plink_glm: thread init failed (phase 2): %s", errstr_buf);
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
		throw IOException("plink_glm: PgrInit failed for '%s'", bind_data.pgen_path);
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

	// Pre-allocate logistic regression scratch buffers (once per thread)
	if (bind_data.is_logistic) {
		state->logistic_bufs.Allocate(bind_data.effective_sample_ct, bind_data.predictor_ct, bind_data.use_firth);
	}

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Per-variant regression results
// ---------------------------------------------------------------------------

struct GlmResult {
	double beta = NAN;
	double se = NAN;
	double t_stat = NAN;
	double p_value = NAN;
	double a1_freq = NAN;
	double odds_ratio = NAN;
	uint32_t obs_ct = 0;
	const char *errcode = nullptr;
	bool firth_applied = false;
	bool is_logistic = false;
};

// ---------------------------------------------------------------------------
// Linear regression (OLS)
//
// Handles both simple (no covariates) and multivariate (with covariates).
// Missing genotypes (-9.0) and missing phenotypes (NaN) are excluded
// per-variant, so OBS_CT may vary across variants.
// ---------------------------------------------------------------------------

static GlmResult ComputeLinearRegression(const double *dosages, const double *phenotype,
                                         const vector<vector<double>> &covariates, uint32_t sample_ct) {
	GlmResult result;

	int n_covars = static_cast<int>(covariates.size());
	int p = 2 + n_covars; // intercept + genotype + covariates

	// First pass: count non-missing and compute allele frequency
	uint32_t n = 0;
	double sum_x = 0.0;
	for (uint32_t i = 0; i < sample_ct; i++) {
		if (dosages[i] == -9.0 || std::isnan(phenotype[i])) {
			continue;
		}
		sum_x += dosages[i];
		n++;
	}

	result.obs_ct = n;

	if (n < static_cast<uint32_t>(p + 1)) {
		result.errcode = "TOO_FEW_SAMPLES";
		return result;
	}

	result.a1_freq = sum_x / (2.0 * static_cast<double>(n));

	// Simple case: no covariates — use optimized formula
	if (n_covars == 0) {
		double sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0, sum_yy = 0.0;
		for (uint32_t i = 0; i < sample_ct; i++) {
			double x = dosages[i];
			double y = phenotype[i];
			if (x == -9.0 || std::isnan(y)) {
				continue;
			}
			sum_y += y;
			sum_xx += x * x;
			sum_xy += x * y;
			sum_yy += y * y;
		}

		double nd = static_cast<double>(n);
		double sxx = sum_xx - sum_x * sum_x / nd;
		double sxy = sum_xy - sum_x * sum_y / nd;
		double syy = sum_yy - sum_y * sum_y / nd;

		if (sxx < 1e-20) {
			result.errcode = "CONST_ALLELE";
			return result;
		}

		result.beta = sxy / sxx;
		double rss = syy - sxy * sxy / sxx;
		if (rss < 0.0) {
			rss = 0.0;
		}

		double df = nd - 2.0;
		double mse = rss / df;
		double se_sq = mse / sxx;
		if (se_sq < 1e-30) {
			result.errcode = "ZERO_VARIANCE";
			return result;
		}

		result.se = std::sqrt(se_sq);
		result.t_stat = result.beta / result.se;
		result.p_value = TstatToPvalue(result.t_stat, df);
		return result;
	}

	// Multivariate case: use plink2's LinearRegressionInv
	uint32_t up = static_cast<uint32_t>(p);

	// Build predictor-major double matrix: predictors_pmaj[p × n]
	// Row 0: intercept (all 1.0)
	// Row 1: genotype dosages (non-missing only)
	// Rows 2+: covariates (non-missing only)
	vector<double> predictors_pmaj(up * n);
	vector<double> pheno_d(n);
	double pheno_ssq = 0.0;

	uint32_t nm_idx = 0;
	for (uint32_t i = 0; i < sample_ct; i++) {
		if (dosages[i] == -9.0 || std::isnan(phenotype[i])) {
			continue;
		}
		predictors_pmaj[0 * n + nm_idx] = 1.0;        // intercept
		predictors_pmaj[1 * n + nm_idx] = dosages[i]; // genotype
		for (int c = 0; c < n_covars; c++) {
			predictors_pmaj[(2 + c) * n + nm_idx] = covariates[c][i];
		}
		pheno_d[nm_idx] = phenotype[i];
		pheno_ssq += phenotype[i] * phenotype[i];
		nm_idx++;
	}

	// Check for constant genotype
	double geno_mean = sum_x / static_cast<double>(n);
	double geno_var = 0.0;
	for (uint32_t idx = 0; idx < n; idx++) {
		double d = predictors_pmaj[1 * n + idx] - geno_mean;
		geno_var += d * d;
	}
	if (geno_var < 1e-20) {
		result.errcode = "CONST_ALLELE";
		return result;
	}

	// Compute X'X (predictors_pmaj * transpose)
	vector<double> xtx_inv(up * up);
	MultiplySelfTranspose(predictors_pmaj.data(), up, n, xtx_inv.data());

	// Allocate workspace
	vector<double> fitted_coefs(up);
	vector<double> xt_y(up);
	vector<MatrixInvertBuf1> mi_buf(2 * up);
	vector<double> dbl_2d_buf(up * std::max(up, 7u));

	// Call plink2's LinearRegressionInv
	if (LinearRegressionInv(pheno_d.data(), predictors_pmaj.data(), up, n, 1, xtx_inv.data(), fitted_coefs.data(),
	                        xt_y.data(), mi_buf.data(), dbl_2d_buf.data())) {
		result.errcode = "SINGULAR_MATRIX";
		return result;
	}

	// RSS = pheno_ssq - dot(xt_y, fitted_coefs)
	double rss = pheno_ssq;
	for (uint32_t j = 0; j < up; j++) {
		rss -= xt_y[j] * fitted_coefs[j];
	}
	if (rss < 0.0) {
		rss = 0.0;
	}

	double df = static_cast<double>(n) - static_cast<double>(p);
	double mse = rss / df;

	// Extract beta and SE for genotype coefficient (index 1)
	// xtx_inv is now (X'X)^{-1}, row-major p×p
	result.beta = fitted_coefs[1];
	double se_sq = mse * xtx_inv[1 * up + 1];
	if (se_sq < 1e-30) {
		result.errcode = "ZERO_VARIANCE";
		return result;
	}
	result.se = std::sqrt(se_sq);
	result.t_stat = result.beta / result.se;
	result.p_value = TstatToPvalue(result.t_stat, df);

	return result;
}

// ---------------------------------------------------------------------------
// Logistic regression (IRLS with optional Firth correction)
// ---------------------------------------------------------------------------

static GlmResult ComputeLogisticRegression(const double *dosages, const double *phenotype,
                                           const vector<vector<double>> &covariates, uint32_t sample_ct, bool use_firth,
                                           LogisticBuffers &bufs) {
	GlmResult result;
	result.is_logistic = true;

	int n_covars = static_cast<int>(covariates.size());
	uint32_t p = static_cast<uint32_t>(2 + n_covars); // intercept + genotype + covariates

	// Collect non-missing samples into pre-allocated buffer
	uint32_t nm_count = 0;
	double sum_x = 0.0;
	for (uint32_t i = 0; i < sample_ct; i++) {
		if (dosages[i] == -9.0 || std::isnan(phenotype[i])) {
			continue;
		}
		bufs.nm_indices[nm_count++] = i;
		sum_x += dosages[i];
	}

	uint32_t n = nm_count;
	result.obs_ct = n;

	if (n < p + 1) {
		result.errcode = "TOO_FEW_SAMPLES";
		return result;
	}

	result.a1_freq = sum_x / (2.0 * static_cast<double>(n));

	// Check for constant genotype
	double geno_mean = sum_x / static_cast<double>(n);
	double geno_var = 0.0;
	for (uint32_t idx = 0; idx < n; idx++) {
		double d = dosages[bufs.nm_indices[idx]] - geno_mean;
		geno_var += d * d;
	}
	if (geno_var < 1e-20) {
		result.errcode = "CONST_ALLELE";
		return result;
	}

	// Vector-aligned sizes for plink2's SIMD code
	uint32_t sample_ctav = RoundUpPow2(n, kFloatPerFVec);
	uint32_t predictor_ctav = RoundUpPow2(p, kFloatPerFVec);

	// Zero scratch buffers for this variant's sample count
	bufs.Reset(n, p, use_firth);

	// Build covariate-major float matrix xx[p × sample_ctav]:
	// Row 0: intercept (1.0f), zero-padded to sample_ctav
	// Row 1: genotype dosages (float), zero-padded
	// Rows 2+: covariates (float), zero-padded
	for (uint32_t idx = 0; idx < n; idx++) {
		uint32_t i = bufs.nm_indices[idx];
		bufs.xx[0 * sample_ctav + idx] = 1.0f;                           // intercept
		bufs.xx[1 * sample_ctav + idx] = static_cast<float>(dosages[i]); // genotype
		for (int c = 0; c < n_covars; c++) {
			bufs.xx[(2 + c) * sample_ctav + idx] = static_cast<float>(covariates[c][i]);
		}
		bufs.yy[idx] = static_cast<float>(phenotype[i]);
	}

	uint32_t is_unfinished = 0;
	plink2::BoolErr logistic_failed = plink2::LogisticRegressionF(
	    bufs.yy.data(), bufs.xx.data(), nullptr, n, p, bufs.coef.data(), &is_unfinished, bufs.ll.data(), bufs.pp.data(),
	    bufs.vv.data(), bufs.hh.data(), bufs.grad.data(), bufs.dcoef.data());

	// Treat "unfinished" (hit iteration limit) as failure for Firth fallback
	bool needs_firth_fallback = logistic_failed || is_unfinished;
	bool firth_applied = false;

	if (needs_firth_fallback && use_firth) {
		// Logistic failed — try Firth correction using pre-allocated buffers
		std::fill_n(bufs.coef.data(), predictor_ctav, 0.0f);
		std::fill_n(bufs.half_inverted_buf.data(), p * p, 0.0);
		std::fill_n(bufs.dbl_2d_buf.data(), p * p, 0.0);
		std::fill_n(bufs.ustar.data(), predictor_ctav, 0.0f);
		std::fill_n(bufs.delta.data(), predictor_ctav, 0.0f);
		std::fill_n(bufs.hdiag.data(), sample_ctav, 0.0f);
		std::fill_n(bufs.ww.data(), sample_ctav, 0.0f);
		std::fill_n(bufs.hh0.data(), p * predictor_ctav, 0.0f);
		std::fill_n(bufs.tmpnxk_buf.data(), p * sample_ctav, 0.0f);
		is_unfinished = 0;

		plink2::BoolErr firth_failed = plink2::FirthRegressionF(
		    bufs.yy.data(), bufs.xx.data(), nullptr, n, p, bufs.coef.data(), &is_unfinished, bufs.hh.data(),
		    bufs.half_inverted_buf.data(), bufs.inv_1d_buf.data(), bufs.dbl_2d_buf.data(), bufs.pp.data(),
		    bufs.vv.data(), bufs.ustar.data(), bufs.delta.data(), bufs.hdiag.data(), bufs.ww.data(), bufs.hh0.data(),
		    bufs.tmpnxk_buf.data());

		if (firth_failed) {
			result.errcode = "NO_CONVERGENCE";
			return result;
		}
		firth_applied = true;
		// After Firth, hh contains the inverted variance-covariance matrix
	} else if (needs_firth_fallback) {
		// Logistic failed/unfinished, no Firth
		result.errcode = logistic_failed ? "SEPARATION" : "NO_CONVERGENCE";
		return result;
	} else {
		// Logistic succeeded (with or without Firth enabled) — need to invert hh for SE.
		// hh contains the Hessian from the last IRLS iteration.
		// Note: we don't call ReflectFmatrix0 after inversion because we only read
		// diagonal elements (hh[i * predictor_ctav + i]). Off-diagonal elements in
		// the upper triangle are not filled — callers must not rely on them.
		std::fill_n(bufs.half_inverted_buf.data(), p * p, 0.0);
		std::fill_n(bufs.dbl_2d_buf.data(), p * p, 0.0);
		if (InvertSymmdefFmatrixFirstHalf(p, predictor_ctav, bufs.hh.data(), bufs.half_inverted_buf.data(),
		                                  bufs.inv_1d_buf.data(), bufs.dbl_2d_buf.data())) {
			result.errcode = "SINGULAR_MATRIX";
			return result;
		}
		InvertSymmdefFmatrixSecondHalf(p, predictor_ctav, bufs.half_inverted_buf.data(), bufs.hh.data(),
		                               bufs.inv_1d_buf.data(), bufs.dbl_2d_buf.data());
	}

	result.firth_applied = firth_applied;

	// Extract results for genotype coefficient (index 1)
	// hh is now inverted variance-covariance, covariate-major with stride predictor_ctav
	result.beta = static_cast<double>(bufs.coef[1]);
	double se_sq = static_cast<double>(bufs.hh[1 * predictor_ctav + 1]);
	if (se_sq < 1e-30) {
		result.errcode = "ZERO_VARIANCE";
		return result;
	}
	result.se = std::sqrt(se_sq);
	result.t_stat = result.beta / result.se; // Wald z-statistic
	result.p_value = ZstatToPvalue(result.t_stat);
	result.odds_ratio = std::exp(result.beta);

	return result;
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static constexpr uint32_t GLM_BATCH_SIZE = 64;

static void PlinkGlmScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkGlmBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkGlmGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkGlmLocalState>();

	auto &column_ids = gstate.column_ids;
	uint32_t end_idx = gstate.end_variant_idx;
	uint32_t sample_ct = bind_data.effective_sample_ct;

	const uintptr_t *sample_include = nullptr;
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		sample_include = bind_data.sample_subset->SampleInclude();
	}

	idx_t rows_emitted = 0;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		uint32_t remaining_capacity = static_cast<uint32_t>(STANDARD_VECTOR_SIZE - rows_emitted);
		uint32_t claim_size = std::min(GLM_BATCH_SIZE, remaining_capacity);
		uint32_t batch_start = gstate.next_variant_idx.fetch_add(claim_size);
		if (batch_start >= end_idx) {
			break;
		}
		uint32_t batch_end = std::min(batch_start + claim_size, end_idx);

		for (uint32_t vidx = batch_start; vidx < batch_end; vidx++) {

			// Compute regression if needed
			GlmResult lr;
			if (gstate.need_regression && lstate.initialized) {
				uint32_t dosage_ct = 0;
				plink2::PglErr err = plink2::PgrGetD(
				    sample_include, lstate.pssi, sample_ct, vidx, &lstate.pgr, lstate.genovec_buf.As<uintptr_t>(),
				    lstate.dosage_present_buf.As<uintptr_t>(), lstate.dosage_main_buf.As<uint16_t>(), &dosage_ct);
				if (err != plink2::kPglRetSuccess) {
					throw IOException("plink_glm: PgrGetD failed for variant %u", vidx);
				}

				plink2::Dosage16ToDoublesMinus9(
				    lstate.genovec_buf.As<uintptr_t>(), lstate.dosage_present_buf.As<uintptr_t>(),
				    lstate.dosage_main_buf.As<uint16_t>(), sample_ct, dosage_ct, lstate.dosage_doubles.data());

				if (bind_data.is_logistic) {
					lr = ComputeLogisticRegression(lstate.dosage_doubles.data(), bind_data.phenotype.data(),
					                               bind_data.covariate_values, sample_ct, bind_data.use_firth,
					                               lstate.logistic_bufs);
				} else {
					lr = ComputeLinearRegression(lstate.dosage_doubles.data(), bind_data.phenotype.data(),
					                             bind_data.covariate_values, sample_ct);
				}
			}

			// Apply p-value threshold filter: skip variants with errors, NaN p-values,
			// or p-values exceeding the threshold
			if (!std::isnan(bind_data.p_threshold)) {
				if (lr.errcode != nullptr || std::isnan(lr.p_value) || lr.p_value > bind_data.p_threshold) {
					continue; // skip this variant, don't emit row
				}
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
					auto val = bind_data.variants.GetChrom(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case COL_POS: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = bind_data.variants.GetPos(vidx);
					break;
				}
				case COL_ID: {
					auto val = bind_data.variants.GetId(vidx);
					if (val.empty()) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_REF: {
					auto val = bind_data.variants.GetRef(vidx);
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					break;
				}
				case COL_ALT: {
					auto val = bind_data.variants.GetAlt(vidx);
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_A1: {
					// Tested allele is ALT (same as plink2)
					auto val = bind_data.variants.GetAlt(vidx);
					if (val.empty() || val == ".") {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, val);
					}
					break;
				}
				case COL_A1_FREQ: {
					if (std::isnan(lr.a1_freq)) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = lr.a1_freq;
					}
					break;
				}
				case COL_TEST: {
					FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, "ADD");
					break;
				}
				case COL_OBS_CT: {
					FlatVector::GetData<int32_t>(vec)[rows_emitted] = static_cast<int32_t>(lr.obs_ct);
					break;
				}
				case COL_BETA: {
					if (lr.errcode != nullptr) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = lr.beta;
					}
					break;
				}
				case COL_SE: {
					if (lr.errcode != nullptr) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = lr.se;
					}
					break;
				}
				case COL_T_STAT: {
					if (lr.errcode != nullptr) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = lr.t_stat;
					}
					break;
				}
				case COL_P: {
					if (lr.errcode != nullptr) {
						FlatVector::SetNull(vec, rows_emitted, true);
					} else {
						FlatVector::GetData<double>(vec)[rows_emitted] = lr.p_value;
					}
					break;
				}
				case COL_ERRCODE: {
					if (lr.errcode != nullptr) {
						FlatVector::GetData<string_t>(vec)[rows_emitted] = StringVector::AddString(vec, lr.errcode);
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
					break;
				}
				case COL_OR: {
					if (lr.is_logistic && lr.errcode == nullptr) {
						FlatVector::GetData<double>(vec)[rows_emitted] = lr.odds_ratio;
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
					}
					break;
				}
				case COL_FIRTH_YN: {
					if (lr.is_logistic && lr.errcode == nullptr) {
						FlatVector::GetData<string_t>(vec)[rows_emitted] =
						    StringVector::AddString(vec, lr.firth_applied ? "Y" : "N");
					} else {
						FlatVector::SetNull(vec, rows_emitted, true);
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

void RegisterPlinkGlm(ExtensionLoader &loader) {
	TableFunction plink_glm("plink_glm", {LogicalType::VARCHAR}, PlinkGlmScan, PlinkGlmBind, PlinkGlmInitGlobal,
	                        PlinkGlmInitLocal);

	plink_glm.projection_pushdown = true;

	plink_glm.named_parameters["pvar"] = LogicalType::VARCHAR;
	plink_glm.named_parameters["psam"] = LogicalType::VARCHAR;
	plink_glm.named_parameters["phenotype"] = LogicalType::ANY;
	plink_glm.named_parameters["covariates"] = LogicalType::ANY;
	plink_glm.named_parameters["samples"] = LogicalType::ANY;
	plink_glm.named_parameters["region"] = LogicalType::VARCHAR;
	plink_glm.named_parameters["model"] = LogicalType::VARCHAR;
	plink_glm.named_parameters["firth"] = LogicalType::BOOLEAN;
	plink_glm.named_parameters["p_threshold"] = LogicalType::DOUBLE;

	loader.RegisterFunction(plink_glm);
}

} // namespace duckdb
