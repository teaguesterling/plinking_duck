#include "plink_ld.hpp"
#include "duckdb_compat.hpp"
#include "plink_common.hpp"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <unordered_map>
#include <vector>

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
static constexpr idx_t COL_R = 6;
static constexpr idx_t COL_R2 = 7;
static constexpr idx_t COL_D_PRIME = 8;
static constexpr idx_t COL_OBS_CT = 9;

// Per-thread genotype caching prevents the windowed LD scan from rereading the
// same partner variants thousands of times. Keep this bounded so whole-genome
// or very large-region scans fall back to streaming PgrGet() calls.
static constexpr uintptr_t LD_CACHE_MAX_BYTES_PER_THREAD = 128ULL * 1024ULL * 1024ULL;
static constexpr uintptr_t LD_CACHE_MAX_BYTES_TOTAL = 512ULL * 1024ULL * 1024ULL;

// ---------------------------------------------------------------------------
// String / population-weight helpers
// ---------------------------------------------------------------------------

static string TrimCopy(const string &value) {
	idx_t start = 0;
	idx_t end = value.size();
	while (start < end && std::isspace(static_cast<unsigned char>(value[start]))) {
		start++;
	}
	while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
		end--;
	}
	return value.substr(start, end - start);
}

static string UpperCopy(string value) {
	std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) { return std::toupper(ch); });
	return value;
}

static vector<string> SplitDelimited(const string &value, char delim) {
	vector<string> out;
	std::stringstream ss(value);
	string item;
	while (std::getline(ss, item, delim)) {
		out.push_back(item);
	}
	return out;
}

struct LdRegionSpec {
	string chrom;
	int64_t start = 0;
	int64_t end = 0;
	bool has_region = false;
};

static LdRegionSpec ParseLdRegionSpec(const string &region_str) {
	auto colon = region_str.find(':');
	if (colon == string::npos || colon == 0) {
		throw InvalidInputException("plink_ld: invalid region format '%s' (expected 'chr:start-end')", region_str);
	}
	string chrom = region_str.substr(0, colon);
	string range_part = region_str.substr(colon + 1);
	auto dash = range_part.find('-');
	if (dash == string::npos) {
		throw InvalidInputException("plink_ld: invalid region format '%s' (expected 'chr:start-end')", region_str);
	}
	auto parse_pos = [&](const string &text, const char *label) -> int64_t {
		char *parse_end = nullptr;
		errno = 0;
		long long value = std::strtoll(text.c_str(), &parse_end, 10);
		if (parse_end == text.c_str() || *parse_end != '\0' || errno != 0 || value < 0) {
			throw InvalidInputException("plink_ld: invalid region %s position in '%s'", label, region_str);
		}
		return static_cast<int64_t>(value);
	};
	LdRegionSpec spec;
	spec.chrom = chrom;
	spec.start = parse_pos(range_part.substr(0, dash), "start");
	spec.end = parse_pos(range_part.substr(dash + 1), "end");
	if (spec.start > spec.end) {
		throw InvalidInputException("plink_ld: region start is greater than end in '%s'", region_str);
	}
	spec.has_region = true;
	return spec;
}

static unordered_map<string, double> ParsePopulationWeightsSpec(const string &spec) {
	unordered_map<string, double> weights;
	double total = 0.0;
	for (auto entry : SplitDelimited(spec, ',')) {
		entry = TrimCopy(entry);
		if (entry.empty()) {
			continue;
		}
		auto eq = entry.find('=');
		if (eq == string::npos) {
			eq = entry.find(':');
		}
		if (eq == string::npos || eq == 0 || eq == entry.size() - 1) {
			throw InvalidInputException("plink_ld: population_weights entries must be 'POP=weight' or 'POP:weight'");
		}
		string pop = UpperCopy(TrimCopy(entry.substr(0, eq)));
		string weight_str = TrimCopy(entry.substr(eq + 1));
		char *parse_end = nullptr;
		errno = 0;
		double weight = std::strtod(weight_str.c_str(), &parse_end);
		if (parse_end == weight_str.c_str() || *parse_end != '\0' || errno != 0 || !std::isfinite(weight) ||
		    weight < 0.0) {
			throw InvalidInputException("plink_ld: invalid non-negative population weight '%s'", weight_str);
		}
		weights[pop] += weight;
		total += weight;
	}
	if (weights.empty() || total <= 0.0) {
		throw InvalidInputException("plink_ld: population_weights must contain at least one positive weight");
	}
	for (auto &kv : weights) {
		kv.second /= total;
	}
	return weights;
}

// ---------------------------------------------------------------------------
// LD computation
// ---------------------------------------------------------------------------

enum class LdMode : uint8_t { PAIRWISE, WINDOWED };

struct LdResult {
	double r;
	double r2;
	double d_prime;
	uint32_t obs_ct;
	bool is_valid; // false if monomorphic, < 2 obs, etc.
};

struct GenovecStats {
	double sum = 0;
	double sum2 = 0;
	uint32_t obs_ct = 0;
	bool no_missing = true;
};

static inline uint32_t PopcountWord(uintptr_t x) {
#if UINTPTR_MAX == UINT64_MAX
	return static_cast<uint32_t>(__builtin_popcountll(static_cast<unsigned long long>(x)));
#else
	return static_cast<uint32_t>(__builtin_popcount(static_cast<unsigned int>(x)));
#endif
}

static GenovecStats ComputeGenovecStats(const uintptr_t *genovec, uint32_t sample_ct) {
	GenovecStats stats;
	uint32_t word_ct = plink2::DivUp(sample_ct, plink2::kBitsPerWordD2);
	for (uint32_t widx = 0; widx < word_ct; widx++) {
		uintptr_t word = genovec[widx];
		uint32_t samples_remaining = sample_ct - widx * plink2::kBitsPerWordD2;
		uint32_t samples_in_word = std::min(samples_remaining, static_cast<uint32_t>(plink2::kBitsPerWordD2));
		for (uint32_t sidx = 0; sidx < samples_in_word; sidx++) {
			uint32_t geno = word & 3;
			word >>= 2;
			if (geno == 3) {
				stats.no_missing = false;
				continue;
			}
			stats.sum += static_cast<double>(geno);
			stats.sum2 += static_cast<double>(geno * geno);
			stats.obs_ct++;
		}
	}
	return stats;
}

static inline uintptr_t MaskTrailingGenotypes(uintptr_t word, uint32_t samples_in_word) {
	if (samples_in_word >= static_cast<uint32_t>(plink2::kBitsPerWordD2)) {
		return word;
	}
	uintptr_t valid_bits = static_cast<uintptr_t>(samples_in_word) * 2;
	uintptr_t valid_mask = (static_cast<uintptr_t>(1) << valid_bits) - 1;
	return word & valid_mask;
}

static double DotNoMissingGenovec(const uintptr_t *genovec_a, const uintptr_t *genovec_b, uint32_t sample_ct) {
	static constexpr uintptr_t even_bit_mask = ~static_cast<uintptr_t>(0) / 3;
	double sum_ab = 0;
	uint32_t word_ct = plink2::DivUp(sample_ct, plink2::kBitsPerWordD2);
	for (uint32_t widx = 0; widx < word_ct; widx++) {
		uintptr_t word_a = genovec_a[widx];
		uintptr_t word_b = genovec_b[widx];
		uint32_t samples_remaining = sample_ct - widx * plink2::kBitsPerWordD2;
		uint32_t samples_in_word = std::min(samples_remaining, static_cast<uint32_t>(plink2::kBitsPerWordD2));
		word_a = MaskTrailingGenotypes(word_a, samples_in_word);
		word_b = MaskTrailingGenotypes(word_b, samples_in_word);

		uintptr_t lo_a = word_a & even_bit_mask;
		uintptr_t hi_a = (word_a >> 1) & even_bit_mask;
		uintptr_t lo_b = word_b & even_bit_mask;
		uintptr_t hi_b = (word_b >> 1) & even_bit_mask;

		sum_ab += static_cast<double>(PopcountWord(lo_a & lo_b));
		sum_ab += 2.0 * static_cast<double>(PopcountWord(lo_a & hi_b));
		sum_ab += 2.0 * static_cast<double>(PopcountWord(hi_a & lo_b));
		sum_ab += 4.0 * static_cast<double>(PopcountWord(hi_a & hi_b));
	}
	return sum_ab;
}

static LdResult FinalizeLdStats(double sum_a, double sum_b, double sum_ab, double sum_a2, double sum_b2, uint32_t n,
                                bool compute_r, bool compute_d_prime) {
	LdResult result;
	result.obs_ct = n;
	result.is_valid = false;
	result.r = 0;
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
	double var_prod = var_a * var_b;
	if (compute_r) {
		result.r = cov_ab / std::sqrt(var_prod);
		result.r2 = result.r * result.r;
	} else {
		result.r2 = (cov_ab * cov_ab) / var_prod;
	}

	if (!compute_d_prime) {
		return result;
	}

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

static LdResult ComputeLdStatsNoMissing(const uintptr_t *genovec_a, const uintptr_t *genovec_b, uint32_t sample_ct,
                                        const GenovecStats &stats_a, const GenovecStats &stats_b, bool compute_r,
                                        bool compute_d_prime) {
	double sum_ab = DotNoMissingGenovec(genovec_a, genovec_b, sample_ct);
	return FinalizeLdStats(stats_a.sum, stats_b.sum, sum_ab, stats_a.sum2, stats_b.sum2, sample_ct, compute_r,
	                       compute_d_prime);
}

//! Compute LD statistics from two packed 2-bit genotype arrays.
//! Genotype encoding: 0=hom_ref, 1=het, 2=hom_alt, 3=missing.
static LdResult ComputeLdStats(const uintptr_t *genovec_a, const uintptr_t *genovec_b, uint32_t sample_ct,
                               bool compute_r, bool compute_d_prime) {
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
	result.r = 0;
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
	double var_prod = var_a * var_b;
	if (compute_r) {
		result.r = cov_ab / std::sqrt(var_prod);
		result.r2 = result.r * result.r;
	} else {
		result.r2 = (cov_ab * cov_ab) / var_prod;
	}

	if (!compute_d_prime) {
		return result;
	}

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

static LdResult ComputeWeightedLdStats(const uintptr_t *genovec_a, const uintptr_t *genovec_b, uint32_t sample_ct,
                                       const vector<double> &sample_weights, bool compute_r, bool compute_d_prime) {
	if (sample_weights.size() != sample_ct) {
		throw InternalException("plink_ld: sample weight count does not match effective sample count");
	}

	double sum_w = 0.0, sum_wx = 0.0, sum_wy = 0.0, sum_wxy = 0.0, sum_wx2 = 0.0, sum_wy2 = 0.0;
	uint32_t obs_ct = 0;
	uint32_t sample_idx = 0;
	uint32_t word_ct = plink2::DivUp(sample_ct, plink2::kBitsPerWordD2);
	for (uint32_t widx = 0; widx < word_ct; widx++) {
		uintptr_t word_a = genovec_a[widx];
		uintptr_t word_b = genovec_b[widx];

		uint32_t samples_remaining = sample_ct - widx * plink2::kBitsPerWordD2;
		uint32_t samples_in_word = std::min(samples_remaining, static_cast<uint32_t>(plink2::kBitsPerWordD2));
		for (uint32_t sidx = 0; sidx < samples_in_word; sidx++, sample_idx++) {
			uint32_t geno_a = word_a & 3;
			uint32_t geno_b = word_b & 3;
			word_a >>= 2;
			word_b >>= 2;

			double weight = sample_weights[sample_idx];
			if (weight <= 0.0 || geno_a == 3 || geno_b == 3) {
				continue;
			}

			double ga = static_cast<double>(geno_a);
			double gb = static_cast<double>(geno_b);
			sum_w += weight;
			sum_wx += weight * ga;
			sum_wy += weight * gb;
			sum_wxy += weight * ga * gb;
			sum_wx2 += weight * ga * ga;
			sum_wy2 += weight * gb * gb;
			obs_ct++;
		}
	}

	LdResult result;
	result.obs_ct = obs_ct;
	result.is_valid = false;
	result.r = 0;
	result.r2 = 0;
	result.d_prime = 0;

	if (obs_ct < 2 || sum_w <= 0.0) {
		return result;
	}

	double mean_a = sum_wx / sum_w;
	double mean_b = sum_wy / sum_w;
	double cov_ab = sum_wxy / sum_w - mean_a * mean_b;
	double var_a = sum_wx2 / sum_w - mean_a * mean_a;
	double var_b = sum_wy2 / sum_w - mean_b * mean_b;
	if (var_a < 1e-15 || var_b < 1e-15) {
		return result;
	}

	result.is_valid = true;
	double var_prod = var_a * var_b;
	if (compute_r) {
		result.r = cov_ab / std::sqrt(var_prod);
		result.r2 = result.r * result.r;
	} else {
		result.r2 = (cov_ab * cov_ab) / var_prod;
	}

	if (!compute_d_prime) {
		return result;
	}

	double D = cov_ab / 4.0;
	double p_a = sum_wx / (2.0 * sum_w);
	double p_b = sum_wy / (2.0 * sum_w);
	double D_max;
	if (D >= 0) {
		D_max = std::min(p_a * (1.0 - p_b), (1.0 - p_a) * p_b);
	} else {
		D_max = std::max(-p_a * p_b, -(1.0 - p_a) * (1.0 - p_b));
	}
	if (std::abs(D_max) < 1e-15) {
		result.d_prime = 0.0;
	} else {
		result.d_prime = D / D_max;
	}

	return result;
}

// ---------------------------------------------------------------------------
// Bind data
// ---------------------------------------------------------------------------

struct PlinkLdBindData : public TableFunctionData {
	PlinkLdBindData() {
		plink2::PreinitPgfi(&shared_pgfi);
	}

	~PlinkLdBindData() override {
		if (shared_pgfi_initialized) {
			plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
			plink2::CleanupPgfi(&shared_pgfi, &cleanup_err);
		}
	}

	string pgen_path;
	string pvar_path;
	string psam_path;

	VariantMetadataIndex variants;
	SampleInfo sample_info;
	bool has_sample_info = false;

	uint32_t raw_variant_ct = 0;
	uint32_t raw_sample_ct = 0;

	// Shared immutable pgenlib file metadata. PgfiInitPhase2 can allocate large
	// per-variant index arrays for whole-genome PGENs, so keep one copy in bind
	// data and let each local PgenReader copy the struct/pointers instead of
	// re-running Phase2 per worker.
	plink2::PgenFileInfo shared_pgfi;
	AlignedBuffer shared_pgfi_alloc_buf;
	uint32_t max_vrec_width = 0;
	uintptr_t pgr_alloc_cacheline_ct = 0;
	bool shared_pgfi_initialized = false;

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
	bool cache_genotypes = true;

	// GAUSS-like population weighting. Each included sample receives
	// population_weight / included_population_sample_count, so each population
	// contributes the requested total weight regardless of sample count.
	string population_column;
	string population_weights_spec;
	bool use_population_weights = false;
	vector<double> sample_weights;
};

static vector<string> SplitPsamDataLine(const string &line) {
	if (line.find('\t') != string::npos) {
		return SplitTabLine(line);
	}
	return SplitWhitespaceLine(line);
}

static vector<string> LoadSamplePopulationLabels(ClientContext &context, const string &psam_path,
                                                 const string &population_column, uint32_t raw_sample_ct) {
	if (psam_path.empty()) {
		throw InvalidInputException("plink_ld: population_column requires a .psam/.fam file");
	}
	if (population_column.empty()) {
		throw InvalidInputException("plink_ld: population_column must not be empty when population_weights is supplied");
	}

	auto lines = ReadFileLines(context, psam_path);
	if (lines.empty()) {
		throw InvalidInputException("plink_ld: sample metadata file '%s' is empty", psam_path);
	}

	vector<string> columns;
	idx_t data_start = 0;
	if (!lines[0].empty() && lines[0][0] == '#') {
		columns = SplitPsamDataLine(lines[0]);
		if (!columns.empty() && !columns[0].empty() && columns[0][0] == '#') {
			columns[0] = columns[0].substr(1);
		}
		data_start = 1;
	} else {
		columns = {"FID", "IID", "PAT", "MAT", "SEX", "PHENO"};
	}

	string requested = UpperCopy(population_column);
	idx_t pop_col_idx = columns.size();
	for (idx_t i = 0; i < columns.size(); i++) {
		if (UpperCopy(columns[i]) == requested) {
			pop_col_idx = i;
			break;
		}
	}
	if (pop_col_idx >= columns.size()) {
		throw InvalidInputException("plink_ld: population_column '%s' not found in '%s'", population_column, psam_path);
	}

	vector<string> labels;
	labels.reserve(raw_sample_ct);
	for (idx_t line_idx = data_start; line_idx < lines.size(); line_idx++) {
		auto line = TrimCopy(lines[line_idx]);
		if (line.empty()) {
			continue;
		}
		auto parts = SplitPsamDataLine(line);
		if (pop_col_idx >= parts.size()) {
			throw InvalidInputException("plink_ld: sample metadata row %llu has no population_column '%s'",
			                            static_cast<unsigned long long>(line_idx + 1), population_column);
		}
		labels.push_back(UpperCopy(parts[pop_col_idx]));
	}
	if (labels.size() != raw_sample_ct) {
		throw InvalidInputException("plink_ld: sample metadata '%s' has %llu rows, but .pgen has %u samples", psam_path,
		                            static_cast<unsigned long long>(labels.size()), raw_sample_ct);
	}
	return labels;
}

static vector<uint32_t> IncludedSampleIndices(const PlinkLdBindData &bind_data) {
	vector<uint32_t> included;
	included.reserve(bind_data.effective_sample_ct);
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		const uintptr_t *sample_include = bind_data.sample_subset->SampleInclude();
		for (uint32_t idx = 0; idx < bind_data.raw_sample_ct; idx++) {
			if (plink2::IsSet(sample_include, idx)) {
				included.push_back(idx);
			}
		}
	} else {
		for (uint32_t idx = 0; idx < bind_data.raw_sample_ct; idx++) {
			included.push_back(idx);
		}
	}
	return included;
}

static vector<double> BuildPopulationSampleWeights(const PlinkLdBindData &bind_data, const vector<string> &labels,
                                                   const unordered_map<string, double> &population_weights) {
	auto included = IncludedSampleIndices(bind_data);
	unordered_map<string, uint32_t> included_pop_counts;
	for (auto raw_idx : included) {
		auto it = population_weights.find(labels[raw_idx]);
		if (it != population_weights.end() && it->second > 0.0) {
			included_pop_counts[labels[raw_idx]]++;
		}
	}

	for (auto &kv : population_weights) {
		if (kv.second > 0.0 && included_pop_counts[kv.first] == 0) {
			throw InvalidInputException("plink_ld: population_weights includes '%s', but no included samples have that %s",
			                            kv.first, bind_data.population_column);
		}
	}

	vector<double> sample_weights;
	sample_weights.reserve(included.size());
	double total_weight = 0.0;
	for (auto raw_idx : included) {
		auto it = population_weights.find(labels[raw_idx]);
		double weight = 0.0;
		if (it != population_weights.end() && it->second > 0.0) {
			weight = it->second / static_cast<double>(included_pop_counts[labels[raw_idx]]);
		}
		sample_weights.push_back(weight);
		total_weight += weight;
	}
	if (sample_weights.empty() || total_weight <= 0.0) {
		throw InvalidInputException("plink_ld: population_weights selected zero effective samples");
	}

	// Normalize after dropping unweighted populations so the weighted moments are
	// stable even if the caller specified extra populations outside the current
	// sample subset.
	for (auto &weight : sample_weights) {
		weight /= total_weight;
	}
	return sample_weights;
}

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
	uint32_t max_threads_config = 0;

	// Projection pushdown
	vector<column_t> column_ids;
	bool need_r = false;
	bool need_d_prime = false;
	bool need_any_output = false;

	// Per-local-thread bounded genotype cache settings
	bool use_genotype_cache = false;
	uintptr_t genovec_word_ct = 0;
	uintptr_t genovec_bytes = 0;

	idx_t MaxThreads() const override {
		if (mode == LdMode::PAIRWISE) {
			return 1;
		}
		uint32_t range = end_variant_idx - start_variant_idx;
		return ApplyMaxThreadsCap(range / 50 + 1, max_threads_config);
	}
};

// ---------------------------------------------------------------------------
// Local state (per-thread)
// ---------------------------------------------------------------------------

struct PlinkLdLocalState : public LocalTableFunctionState {
	plink2::PgenReader pgr;
	AlignedBuffer pgr_alloc_buf;

	plink2::PgrSampleSubsetIndex pssi;

	AlignedBuffer genovec_a_buf; // anchor variant genotypes
	AlignedBuffer genovec_b_buf; // partner variant genotypes

	// Optional per-thread dense cache for all genotypes in the active region.
	// This trades bounded memory for eliminating repeated PgrGet() calls in
	// dense window scans.
	AlignedBuffer genotype_cache_buf;
	vector<GenovecStats> genotype_cache_stats;
	bool genotype_cache_ready = false;
	uint32_t cache_start_variant_idx = 0;
	uint32_t cache_end_variant_idx = 0;
	uintptr_t cache_genovec_word_ct = 0;

	// Windowed mode: state preservation across scan calls
	bool in_window = false;
	uint32_t anchor_idx = 0;
	uint32_t next_j = 0;

	bool initialized = false;

	~PlinkLdLocalState() {
		if (initialized) {
			plink2::PglErr reterr = plink2::kPglRetSuccess;
			plink2::CleanupPgr(&pgr, &reterr);
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
		} else if (kv.first == "cache_genotypes") {
			bind_data->cache_genotypes = kv.second.GetValue<bool>();
		} else if (kv.first == "population_column") {
			bind_data->population_column = kv.second.GetValue<string>();
		} else if (kv.first == "population_weights") {
			bind_data->population_weights_spec = kv.second.GetValue<string>();
		} else if (kv.first == "samples" || kv.first == "region") {
			// Handled after pgenlib init
		}
	}

	bind_data->use_population_weights = !bind_data->population_column.empty() || !bind_data->population_weights_spec.empty();
	if (bind_data->use_population_weights &&
	    (bind_data->population_column.empty() || bind_data->population_weights_spec.empty())) {
		throw InvalidInputException("plink_ld: population_column and population_weights must be supplied together");
	}

	LdRegionSpec region_spec;
	auto region_it = input.named_parameters.find("region");
	if (region_it != input.named_parameters.end()) {
		region_spec = ParseLdRegionSpec(region_it->second.GetValue<string>());
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
		bind_data->pvar_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".pvar", ".bim"});
		if (bind_data->pvar_path.empty()) {
			throw InvalidInputException("plink_ld: cannot find .pvar or .bim companion for '%s' "
			                            "(use pvar := 'path' to specify explicitly)",
			                            bind_data->pgen_path);
		}
	}

	if (bind_data->psam_path.empty()) {
		bind_data->psam_path = FindCompanionFileWithParquet(context, fs, bind_data->pgen_path, {".psam", ".fam"});
		// .psam is optional for plink_ld — only needed if samples parameter uses VARCHAR IDs
	}

	// --- Initialize shared pgenlib metadata (Phase 1/2) to get counts and a
	// single immutable per-variant index for all local readers. Re-running Phase
	// 2 per worker duplicates large var_fpos/vrtype arrays on whole-genome PGENs.
	char errstr_buf[plink2::kPglErrstrBufBlen];
	plink2::PgenHeaderCtrl header_ctrl;
	uintptr_t pgfi_alloc_cacheline_ct = 0;

	plink2::PglErr err = plink2::PgfiInitPhase1(bind_data->pgen_path.c_str(), nullptr, UINT32_MAX, UINT32_MAX,
	                                            &header_ctrl, &bind_data->shared_pgfi,
	                                            &pgfi_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&bind_data->shared_pgfi, &cleanup_err);
		throw IOException("plink_ld: failed to open '%s': %s", bind_data->pgen_path, errstr_buf);
	}

	bind_data->raw_variant_ct = bind_data->shared_pgfi.raw_variant_ct;
	bind_data->raw_sample_ct = bind_data->shared_pgfi.raw_sample_ct;

	if (pgfi_alloc_cacheline_ct > 0) {
		bind_data->shared_pgfi_alloc_buf.Allocate(pgfi_alloc_cacheline_ct * plink2::kCacheline);
	}

	err = plink2::PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, bind_data->shared_pgfi.raw_variant_ct,
	                             &bind_data->max_vrec_width, &bind_data->shared_pgfi,
	                             bind_data->shared_pgfi_alloc_buf.As<unsigned char>(),
	                             &bind_data->pgr_alloc_cacheline_ct, errstr_buf);

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgfi(&bind_data->shared_pgfi, &cleanup_err);
		throw IOException("plink_ld: failed to initialize '%s' (phase 2): %s", bind_data->pgen_path, errstr_buf);
	}
	bind_data->shared_pgfi_initialized = true;

	// Force local readers to open their own file handles instead of racing to
	// steal the shared Phase-1 FILE* from the bind-time PgenFileInfo.
	if (bind_data->shared_pgfi.shared_ff) {
		std::fclose(bind_data->shared_pgfi.shared_ff);
		bind_data->shared_pgfi.shared_ff = nullptr;
	}
	if (bind_data->shared_pgfi.pgi_ff) {
		std::fclose(bind_data->shared_pgfi.pgi_ff);
		bind_data->shared_pgfi.pgi_ff = nullptr;
	}

	// --- Load variant metadata ---
	bool metadata_region_pushed = false;
	bool can_region_push_metadata = region_spec.has_region && bind_data->mode == LdMode::WINDOWED;
	if (can_region_push_metadata && IsParquetFile(bind_data->pvar_path)) {
		bind_data->variants = LoadVariantMetadataFromParquetRegion(context, bind_data->pvar_path, region_spec.chrom,
		                                                        region_spec.start, region_spec.end,
		                                                        bind_data->raw_variant_ct, "plink_ld");
		metadata_region_pushed = true;
	} else if (can_region_push_metadata && IsNativePlinkFormat(bind_data->pvar_path)) {
		bind_data->variants = LoadVariantMetadataFromTextRegion(context, bind_data->pvar_path, region_spec.chrom,
		                                                     region_spec.start, region_spec.end, bind_data->raw_variant_ct,
		                                                     "plink_ld");
		metadata_region_pushed = true;
	} else {
		bind_data->variants = LoadVariantMetadata(context, bind_data->pvar_path, "plink_ld");
	}

	if (bind_data->variants.variant_ct != bind_data->raw_variant_ct) {
		throw InvalidInputException("plink_ld: variant count mismatch: .pgen has %u variants, "
		                            ".pvar/.bim '%s' has %llu variants",
		                            bind_data->raw_variant_ct, bind_data->pvar_path,
		                            static_cast<unsigned long long>(bind_data->variants.variant_ct));
	}

	// --- Load sample info (optional) ---
	if (!bind_data->psam_path.empty()) {
		bind_data->sample_info = LoadSampleMetadata(context, bind_data->psam_path);
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

	// --- Process GAUSS-like population weights ---
	if (bind_data->use_population_weights) {
		auto population_weights = ParsePopulationWeightsSpec(bind_data->population_weights_spec);
		auto labels = LoadSamplePopulationLabels(context, bind_data->psam_path, bind_data->population_column,
		                                      bind_data->raw_sample_ct);
		bind_data->sample_weights = BuildPopulationSampleWeights(*bind_data, labels, population_weights);
		if (bind_data->sample_weights.size() != bind_data->effective_sample_ct) {
			throw InternalException("plink_ld: constructed sample weights do not match effective sample count");
		}
	}

	// --- Process region parameter ---
	if (region_spec.has_region) {
		if (metadata_region_pushed) {
			bind_data->variant_range.has_filter = true;
			if (!bind_data->variants.local_to_vidx.empty()) {
				bind_data->variant_range.start_idx = bind_data->variants.local_to_vidx.front();
				bind_data->variant_range.end_idx = bind_data->variants.local_to_vidx.back() + 1;
			}
		} else {
			bind_data->variant_range = ParseRegion(region_it->second.GetValue<string>(), bind_data->variants, "plink_ld");
		}
	}

	// --- Resolve pairwise variant indices ---
	if (bind_data->mode == LdMode::PAIRWISE) {
		// BuildVariantIdIndex handles both dense and sparse (pushdown) indexes,
		// and is O(loaded subset size) — a single pass that's also faster than
		// the previous double-ID search loop.
		auto id_to_idx = BuildVariantIdIndex(bind_data->variants);
		auto it_a = id_to_idx.find(variant1_id);
		if (it_a == id_to_idx.end()) {
			throw InvalidInputException("plink_ld: variant '%s' not found in .pvar", variant1_id);
		}
		bind_data->pairwise_vidx_a = it_a->second;
		auto it_b = id_to_idx.find(variant2_id);
		if (it_b == id_to_idx.end()) {
			throw InvalidInputException("plink_ld: variant '%s' not found in .pvar", variant2_id);
		}
		bind_data->pairwise_vidx_b = it_b->second;
	}

	// --- Register output columns ---
	// R is signed ALT-dosage correlation and is needed by fine-mapping / coloc
	// methods such as SuSiE and ColocBoost. R2 is retained for compatibility.
	names = {"CHROM_A", "POS_A", "ID_A", "CHROM_B", "POS_B", "ID_B", "R", "R2", "D_PRIME", "OBS_CT"};
	return_types = {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::VARCHAR,
	                LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE, LogicalType::INTEGER};

	return std::move(bind_data);
}

// ---------------------------------------------------------------------------
// Init global
// ---------------------------------------------------------------------------

static unique_ptr<GlobalTableFunctionState> PlinkLdInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PlinkLdBindData>();
	auto state = make_uniq<PlinkLdGlobalState>();

	state->mode = bind_data.mode;
	state->max_threads_config = GetPlinkingMaxThreads(context);

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

	state->column_ids = input.column_ids;
	state->need_any_output = false;
	state->need_r = false;
	state->need_d_prime = false;
	for (auto col_id : input.column_ids) {
		if (col_id == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}
		state->need_any_output = true;
		if (col_id == COL_R) {
			state->need_r = true;
		} else if (col_id == COL_D_PRIME) {
			state->need_d_prime = true;
		}
	}

	uint32_t range = state->end_variant_idx - state->start_variant_idx;
	state->genovec_word_ct = plink2::NypCtToAlignedWordCt(bind_data.effective_sample_ct);
	state->genovec_bytes = state->genovec_word_ct * sizeof(uintptr_t);
	uintptr_t cache_bytes = static_cast<uintptr_t>(range) * state->genovec_bytes;
	idx_t estimated_threads = state->mode == LdMode::PAIRWISE ? 1 : ApplyMaxThreadsCap(range / 50 + 1, state->max_threads_config);
	bool cache_fits_per_thread = cache_bytes <= LD_CACHE_MAX_BYTES_PER_THREAD;
	bool cache_fits_total = estimated_threads > 0 && cache_bytes <= (LD_CACHE_MAX_BYTES_TOTAL / estimated_threads);
	state->use_genotype_cache = bind_data.cache_genotypes && state->mode == LdMode::WINDOWED && range > 0 &&
	                            state->genovec_bytes > 0 && cache_fits_per_thread && cache_fits_total;

	return std::move(state);
}

// ---------------------------------------------------------------------------
// Init local (per-thread PgenReader)
// ---------------------------------------------------------------------------

static unique_ptr<LocalTableFunctionState> PlinkLdInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                            GlobalTableFunctionState *global_state) {
	auto &bind_data = input.bind_data->Cast<PlinkLdBindData>();
	auto &gstate = global_state->Cast<PlinkLdGlobalState>();
	auto state = make_uniq<PlinkLdLocalState>();

	// --- Initialize per-thread PgenReader from bind-time shared PgenFileInfo ---
	plink2::PreinitPgr(&state->pgr);
	if (!bind_data.shared_pgfi_initialized) {
		throw InternalException("plink_ld: shared pgen metadata was not initialized");
	}

	if (bind_data.pgr_alloc_cacheline_ct > 0) {
		state->pgr_alloc_buf.Allocate(bind_data.pgr_alloc_cacheline_ct * plink2::kCacheline);
	}

	plink2::PglErr err = plink2::PgrInit(bind_data.pgen_path.c_str(), bind_data.max_vrec_width,
	                                    const_cast<plink2::PgenFileInfo *>(&bind_data.shared_pgfi),
	                                    &state->pgr, state->pgr_alloc_buf.As<unsigned char>());

	if (err != plink2::kPglRetSuccess) {
		plink2::PglErr cleanup_err = plink2::kPglRetSuccess;
		plink2::CleanupPgr(&state->pgr, &cleanup_err);
		throw IOException("plink_ld: PgrInit failed for '%s'", bind_data.pgen_path);
	}

	// Set up sample subsetting
	if (bind_data.has_sample_subset && bind_data.sample_subset) {
		plink2::PgrSetSampleSubsetIndex(bind_data.sample_subset->CumulativePopcounts(), &state->pgr, &state->pssi);
	} else {
		plink2::PgrClearSampleSubsetIndex(&state->pgr, &state->pssi);
	}

	// Allocate two genovec buffers (anchor + partner)
	uintptr_t genovec_word_ct = gstate.genovec_word_ct ? gstate.genovec_word_ct : plink2::NypCtToAlignedWordCt(bind_data.effective_sample_ct);
	uintptr_t genovec_bytes = gstate.genovec_bytes ? gstate.genovec_bytes : genovec_word_ct * sizeof(uintptr_t);
	state->genovec_a_buf.Allocate(genovec_bytes);
	state->genovec_b_buf.Allocate(genovec_bytes);
	std::memset(state->genovec_a_buf.ptr, 0, genovec_bytes);
	std::memset(state->genovec_b_buf.ptr, 0, genovec_bytes);

	if (gstate.use_genotype_cache) {
		uint32_t range = gstate.end_variant_idx - gstate.start_variant_idx;
		uintptr_t cache_bytes = static_cast<uintptr_t>(range) * genovec_bytes;
		state->genotype_cache_buf.Allocate(cache_bytes);
		state->cache_start_variant_idx = gstate.start_variant_idx;
		state->cache_end_variant_idx = gstate.end_variant_idx;
		state->cache_genovec_word_ct = genovec_word_ct;
		auto *cache = state->genotype_cache_buf.As<uintptr_t>();
		state->genotype_cache_stats.reserve(range);
		for (uint32_t vidx = state->cache_start_variant_idx; vidx < state->cache_end_variant_idx; vidx++) {
			auto *dest = cache + static_cast<uintptr_t>(vidx - state->cache_start_variant_idx) * genovec_word_ct;
			const uintptr_t *sample_include = nullptr;
			if (bind_data.has_sample_subset && bind_data.sample_subset) {
				sample_include = bind_data.sample_subset->SampleInclude();
			}
			plink2::PglErr cache_err = plink2::PgrGet(sample_include, state->pssi, bind_data.effective_sample_ct, vidx,
			                                         &state->pgr, dest);
			if (cache_err != plink2::kPglRetSuccess) {
				throw IOException("plink_ld: PgrGet failed while caching variant %u", vidx);
			}
			state->genotype_cache_stats.push_back(ComputeGenovecStats(dest, bind_data.effective_sample_ct));
		}
		state->genotype_cache_ready = true;
	}

	state->initialized = true;
	return std::move(state);
}

// ---------------------------------------------------------------------------
// Helpers: emit an output row and read genotypes
// ---------------------------------------------------------------------------

static void EmitRow(DataChunk &output, idx_t row_idx, const PlinkLdBindData &bind_data,
                    const PlinkLdGlobalState &gstate, uint32_t vidx_a, uint32_t vidx_b, const LdResult &result) {
	auto &variants = bind_data.variants;

	for (idx_t out_col = 0; out_col < gstate.column_ids.size(); out_col++) {
		auto file_col = gstate.column_ids[out_col];
		if (file_col == COLUMN_IDENTIFIER_ROW_ID) {
			continue;
		}
		auto &vec = output.data[out_col];
		switch (file_col) {
		case COL_CHROM_A:
			FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, variants.GetChrom(vidx_a));
			break;
		case COL_POS_A:
			FlatVector::GetData<int32_t>(vec)[row_idx] = variants.GetPos(vidx_a);
			break;
		case COL_ID_A: {
			auto id_a = variants.GetId(vidx_a);
			if (id_a.empty()) {
				FlatVector::SetNull(vec, row_idx, true);
			} else {
				FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, id_a);
			}
			break;
		}
		case COL_CHROM_B:
			FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, variants.GetChrom(vidx_b));
			break;
		case COL_POS_B:
			FlatVector::GetData<int32_t>(vec)[row_idx] = variants.GetPos(vidx_b);
			break;
		case COL_ID_B: {
			auto id_b = variants.GetId(vidx_b);
			if (id_b.empty()) {
				FlatVector::SetNull(vec, row_idx, true);
			} else {
				FlatVector::GetData<string_t>(vec)[row_idx] = StringVector::AddString(vec, id_b);
			}
			break;
		}
		case COL_R:
			if (result.is_valid) {
				FlatVector::GetData<double>(vec)[row_idx] = result.r;
			} else {
				FlatVector::SetNull(vec, row_idx, true);
			}
			break;
		case COL_R2:
			if (result.is_valid) {
				FlatVector::GetData<double>(vec)[row_idx] = result.r2;
			} else {
				FlatVector::SetNull(vec, row_idx, true);
			}
			break;
		case COL_D_PRIME:
			if (result.is_valid) {
				FlatVector::GetData<double>(vec)[row_idx] = result.d_prime;
			} else {
				FlatVector::SetNull(vec, row_idx, true);
			}
			break;
		case COL_OBS_CT:
			FlatVector::GetData<int32_t>(vec)[row_idx] = static_cast<int32_t>(result.obs_ct);
			break;
		default:
			break;
		}
	}
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

static const uintptr_t *GetGenovec(PlinkLdLocalState &lstate, const PlinkLdBindData &bind_data, uint32_t vidx,
                                   uintptr_t *fallback_buf) {
	if (lstate.genotype_cache_ready && vidx >= lstate.cache_start_variant_idx && vidx < lstate.cache_end_variant_idx) {
		auto *cache = lstate.genotype_cache_buf.As<uintptr_t>();
		return cache + static_cast<uintptr_t>(vidx - lstate.cache_start_variant_idx) * lstate.cache_genovec_word_ct;
	}
	ReadGenovec(lstate, bind_data, vidx, fallback_buf);
	return fallback_buf;
}

static const GenovecStats *GetCachedGenovecStats(const PlinkLdLocalState &lstate, uint32_t vidx) {
	if (lstate.genotype_cache_ready && vidx >= lstate.cache_start_variant_idx && vidx < lstate.cache_end_variant_idx) {
		return &lstate.genotype_cache_stats[vidx - lstate.cache_start_variant_idx];
	}
	return nullptr;
}

static LdResult ComputeLdStatsDispatch(const PlinkLdBindData &bind_data, const uintptr_t *genovec_a,
                                       const uintptr_t *genovec_b, uint32_t sample_ct, const GenovecStats *stats_a,
                                       const GenovecStats *stats_b, bool compute_r, bool compute_d_prime) {
	if (bind_data.use_population_weights) {
		return ComputeWeightedLdStats(genovec_a, genovec_b, sample_ct, bind_data.sample_weights, compute_r,
		                              compute_d_prime);
	}
	if (stats_a && stats_b && stats_a->no_missing && stats_b->no_missing) {
		return ComputeLdStatsNoMissing(genovec_a, genovec_b, sample_ct, *stats_a, *stats_b, compute_r, compute_d_prime);
	}
	return ComputeLdStats(genovec_a, genovec_b, sample_ct, compute_r, compute_d_prime);
}

// ---------------------------------------------------------------------------
// Scan function
// ---------------------------------------------------------------------------

static void PlinkLdScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<PlinkLdBindData>();
	auto &gstate = data_p.global_state->Cast<PlinkLdGlobalState>();
	auto &lstate = data_p.local_state->Cast<PlinkLdLocalState>();

	if (!lstate.initialized) {
		CompatSetOutputCardinality(output, 0);
		return;
	}

	uint32_t sample_ct = bind_data.effective_sample_ct;
	auto *genovec_a_buf = lstate.genovec_a_buf.As<uintptr_t>();
	auto *genovec_b_buf = lstate.genovec_b_buf.As<uintptr_t>();
	const uintptr_t *genovec_a = nullptr;
	const uintptr_t *genovec_b = nullptr;
	const GenovecStats *stats_a = nullptr;
	const GenovecStats *stats_b = nullptr;

	idx_t rows_emitted = 0;

	if (gstate.mode == LdMode::PAIRWISE) {
		// --- Pairwise mode: emit a single pair ---
		if (gstate.pair_emitted.exchange(true)) {
			CompatSetOutputCardinality(output, 0);
			return;
		}

		uint32_t vidx_a = bind_data.pairwise_vidx_a;
		uint32_t vidx_b = bind_data.pairwise_vidx_b;

		genovec_a = GetGenovec(lstate, bind_data, vidx_a, genovec_a_buf);
		stats_a = GetCachedGenovecStats(lstate, vidx_a);
		if (vidx_a == vidx_b) {
			// Self-LD: use same buffer for both
			auto result = ComputeLdStatsDispatch(bind_data, genovec_a, genovec_a, sample_ct, stats_a, stats_a,
			                                    gstate.need_r, gstate.need_d_prime);
			EmitRow(output, 0, bind_data, gstate, vidx_a, vidx_b, result);
		} else {
			genovec_b = GetGenovec(lstate, bind_data, vidx_b, genovec_b_buf);
			stats_b = GetCachedGenovecStats(lstate, vidx_b);
			auto result = ComputeLdStatsDispatch(bind_data, genovec_a, genovec_b, sample_ct, stats_a, stats_b,
			                                    gstate.need_r, gstate.need_d_prime);
			EmitRow(output, 0, bind_data, gstate, vidx_a, vidx_b, result);
		}
		CompatSetOutputCardinality(output, 1);
		return;
	}

	// --- Windowed mode ---
	uint32_t end_idx = gstate.end_variant_idx;
	auto &variants = bind_data.variants;

	// Cache anchor chrom/pos for the inner loop to avoid repeated parsing
	string anchor_chrom;
	int32_t anchor_pos = 0;
	bool anchor_cached = false;

	while (rows_emitted < STANDARD_VECTOR_SIZE) {
		if (lstate.in_window) {
			// Resume scanning partners for current anchor
			uint32_t ai = lstate.anchor_idx;
			uint32_t j = lstate.next_j;
			genovec_a = GetGenovec(lstate, bind_data, ai, genovec_a_buf);
			stats_a = GetCachedGenovecStats(lstate, ai);

			if (!anchor_cached) {
				anchor_chrom = variants.GetChrom(ai);
				anchor_pos = variants.GetPos(ai);
				anchor_cached = true;
			}

			while (j < end_idx) {
				auto j_chrom = variants.GetChrom(j);
				bool same_chrom = (j_chrom == anchor_chrom);

				if (same_chrom) {
					int64_t dist = static_cast<int64_t>(variants.GetPos(j)) - static_cast<int64_t>(anchor_pos);
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

				genovec_b = GetGenovec(lstate, bind_data, j, genovec_b_buf);
				stats_b = GetCachedGenovecStats(lstate, j);
				auto result = ComputeLdStatsDispatch(bind_data, genovec_a, genovec_b, sample_ct, stats_a, stats_b,
				                                    gstate.need_r, gstate.need_d_prime);

				if (result.is_valid && result.r2 >= bind_data.r2_threshold) {
					EmitRow(output, rows_emitted, bind_data, gstate, ai, j, result);
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
		genovec_a = GetGenovec(lstate, bind_data, anchor_idx, genovec_a_buf);
		stats_a = GetCachedGenovecStats(lstate, anchor_idx);

		lstate.anchor_idx = anchor_idx;
		lstate.next_j = anchor_idx + 1;
		lstate.in_window = true;
		anchor_chrom = variants.GetChrom(anchor_idx);
		anchor_pos = variants.GetPos(anchor_idx);
		anchor_cached = true;
	}

done:
	CompatSetOutputCardinality(output, rows_emitted);
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
	plink_ld.named_parameters["cache_genotypes"] = LogicalType::BOOLEAN;
	plink_ld.named_parameters["population_column"] = LogicalType::VARCHAR;
	plink_ld.named_parameters["population_weights"] = LogicalType::VARCHAR;
	plink_ld.projection_pushdown = true;

	loader.RegisterFunction(plink_ld);
}

} // namespace duckdb
