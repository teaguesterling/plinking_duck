# plink2 / pgenlib Wrapping Audit

**Scope:** Read-only audit of every `src/*.cpp` in the `plinking_duck` DuckDB extension,
scoring where we **wrap** upstream plink2/pgenlib callables vs. where we **reimplement**
their algorithms in our own C++.

**Principle under test:** wrap the plink2 / plink-ng (pgenlib) codebase as much as possible;
reimplementing plink2 logic is a smell to be justified or fixed.

**Method:** For each feature, identified the pgenlib/plink2 calls actually made (evidence of
wrapping), the math actually hand-rolled, and — for each reimplementation — whether a
*callable* upstream equivalent exists and whether it is cleanly linkable or entangled with
plink2's `g_bigstack` arena allocator / `ThreadGroup` / CLI model. That entanglement axis is
the primary discriminator between RED (should have wrapped) and YELLOW (reimplementation forced).

> Note on framing: this is a static, read-only audit. Where I say an upstream file "appears
> cleanly linkable," that means it depends only on already-linked units and shows no
> `bigstack`/thread/CLI references — not that I built it to confirm. Statements are calibrated
> accordingly.

---

## Executive Summary

**Counts:** 🟢 GREEN ×11 · 🟡 YELLOW ×4 · 🔴 RED ×2.

Counting basis: one score per source file (16 `.cpp` files, listed below). `plink_glm.cpp`
is scored **GREEN as a file** (its regression engine wraps plink2) with an embedded **RED
sub-finding** for its p-value layer; that RED is counted separately in the RED tally because it
is the same root cause as the Hardy RED and is the actionable item. So the two REDs are the
Hardy file and the GLM p-value layer; the 11 GREEN + 4 YELLOW cover the 15 remaining files, with
`plink_glm.cpp` appearing in both the GREEN file count and the RED sub-finding. Header files
(`src/include/*.hpp`) are pure declarations — the only inline logic is a trivial 3-element
genotype-mask populator (`plink_common.hpp:451`), no reimplemented statistics — so they are not
scored separately.

The extension is, on the whole, a **faithful wrapper**: all genotype/dosage decoding goes
through pgenlib (`PgrGet`, `PgrGetD`, `PgrGetCounts`, `PgrGetMissingness`,
`GenoarrToBytesMinus9`, `Dosage16ToDoublesMinus9`), the GWAS logistic/Firth solver is a
**verbatim vendored extract** of plink2's own code, and the GLM linear-algebra path **links
and calls** plink2's `plink2_matrix.cc`. The reimplementations that remain fall into two
clean buckets.

### Top items ranked by "most worth moving to a wrap"

1. **🔴 #1 — Re-derived statistics because `plink2_stats.cc` is excluded from the build
   (`CMakeLists.txt:128`).** This is a *single root cause* with two visible symptoms:
   - `plink_hardy.cpp:39 HweExactTest` re-implements the Wigginton HWE exact test, for which
     upstream exposes the callable non-static `HweLnP` (`plink2_stats.cc:1784`) and the
     chrX-specific `HweXchrLnP`.
   - `plink_glm.cpp:136–231` re-implements `BetaIncomplete` / `TstatToPvalue` / `NormalCdf` /
     `ZstatToPvalue`, for which upstream exposes callable `TstatToP2`, `ChisqToP`,
     `QuantileToZscore` (`plink2_stats.h:30,43,51`).

   `plink2_stats.cc` `#include`s only `plink2_stats.h` + `plink2_string.h` (already linked) +
   libc, and a grep for `bigstack|threadpool|ThreadGroup|arena` over it returns **0 hits**. It
   therefore *appears cleanly linkable*. The "intentionally excluded" comment in CMake is the
   smell being audited, not a justification. **One build-list change plausibly unlocks both
   wraps** — highest leverage in the report.

2. **🟡 #2 — `plink2_glm_logistic_math.cpp` (vendored verbatim extract).** A reformatted copy
   of `plink2_glm_logistic.cc` L147–1012 (bodies byte-identical modulo clang-format). Works and
   is faithful, but is a manual-drift burden that must be re-synced whenever the submodule pin
   moves. Extraction was genuinely necessary (home TU is 7106 lines with 77 bigstack/ThreadGroup
   refs). Keep, but track upstream explicitly.

3. **🟡 #3 — `vcf_genotype_parse.cpp` (vendored extract).** Copies
   `VcfConvertUnphasedBiallelicLine` / `VcfConvertPhasedBiallelicLine` / `VcfCheckQuals` out of
   `plink2_import.cc`. Same drift burden; extraction defensible (import.cc is CLI/pipeline-bound).

4. **🟡 #4 — `plink_ld.cpp` genotype-covariance r² / Weir D'.** Forced reimplementation
   (`plink2_ld.cc` is bigstack/CLI-entangled and not linkable), **but** it uses a genotype-level
   Pearson/composite estimator that is numerically *different* from plink2's phased-haplotype r²,
   and its own comments note D' can exceed 1 under HWE deviation. Worth a correctness review even
   though a clean wrap isn't available.

5. **🟡 #5 — `plink_pca.cpp` randomized SVD.** Forced reimplementation (`CalcPca` in
   `plink2_matrix_calc.cc` is bigstack/`PgenMtLoad`/thread-group entangled and not linked).
   Substantial independent math (Eigen `BDCSVD`), so it fails GREEN's "small + faithful" bar,
   but RED would be unfair — no callable equivalent exists.

---

## Per-file findings

### 🔴 `plink_hardy.cpp` — HWE exact test

- **Feature:** per-variant Hardy-Weinberg equilibrium exact-test p-values (ploidy/sex-aware for chrX/Y/MT).
- **Wraps:** genotype counts via `PgrGetCounts` and `PgrGet` (fast, no decompression); standard `Pgfi/Pgr` open/init/cleanup. Sex-aware strata built by counting decoded bytes in `plink_common.cpp:~2090`.
- **Reimplements:** `HweExactTest` (`plink_hardy.cpp:39–129`) — an *independent* implementation of the Wigginton et al. 2005 recurrence (allocates a `vector<double> het_probs`, sums the tail directly).
- **Callable upstream equivalent:** **Yes.** `double HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp)` at `plink2_stats.cc:1784` (non-static; header at `plink2_stats.h:67`) — the same test, by C. Chang, with documented improvements our copy lacks: >64k-genotype overflow safety, no malloc, FP under/overflow handling, log-p output. For chrX, `HweXchrLnP` (`plink2_stats.h:90`) is the proper males-included test; our code instead uses a **female-only** stratum with the autosomal routine (`plink_hardy.cpp:509–535`), a real semantic divergence from plink2.
- **Callability:** `plink2_stats.cc` depends only on `plink2_stats.h` + `plink2_string.h` (linked) + libc; **0** bigstack/thread refs → appears cleanly linkable. It is excluded only by the deliberate `CMakeLists.txt:128` decision.
- **Score:** 🔴 **RED.** A callable, dependency-light, *superior* upstream function exists; the reimplementation is independent (not a faithful copy) and carries overflow + chrX-semantics drift risk. Refactor toward `exp(HweLnP(...))` / `HweXchrLnP(...)` once `plink2_stats.cc` is linked.

### 🔴 `plink_glm.cpp` — p-value / distribution layer

- **Feature:** t- and z-statistic → p-value conversions used to report GWAS `P`.
- **Reimplements:** `LogGamma` (Lanczos, `:120`), `BetaIncomplete` (Lentz continued fraction, `:136`), `TstatToPvalue` (`:202`), `NormalCdf`/`ZstatToPvalue` (`:218–231`). All in `namespace duckdb`, no extraction provenance — original textbook numerics.
- **Callable upstream equivalent:** **Yes.** `TstatToP2` (`plink2_stats.h:43`), `TstatToLnP` (`:47`), `ChisqToP` (`:30`), `QuantileToZscore` (`:51`) — same excluded `plink2_stats.cc`.
- **Score:** 🔴 **RED**, sharing root cause #1. Lower individual stakes than Hardy (standard-normal / t-CDF numerics are low-risk), but it is exactly the "hand-rolled a computation upstream exposes a callable for" pattern. Note: the rest of `plink_glm.cpp` is strongly GREEN (see below) — this is a sub-finding of an otherwise-wrapping file.

### 🟡 `plink2_glm_logistic_math.cpp` — vendored logistic/Firth solver

- **Feature:** float-precision IRLS logistic (`LogisticRegressionF`) + penalized Firth (`FirthRegressionF`) regression and all supporting SIMD kernels; lives in `namespace plink2`.
- **Provenance:** **verbatim extract** of `plink2_glm_logistic.cc` L147–1012. Signatures match upstream exactly (`LogisticRegressionF` upstream L589, `FirthRegressionF` L860, plus `ComputeHessianF`, `CholeskyDecompositionF`, `SolveLinearSystemF`, `LogisticSseF`, etc.); whitespace-normalized bodies show **zero** code differences — deltas are only comment banners and clang-format re-wrapping. Its dense linear algebra is delegated to the linked `plink2_matrix.cc` (`ColMajorFmatrix*`, `InvertSymmdefFmatrix*`).
- **Was extraction necessary?** The two entry points are non-static (callable in principle), but their home TU is 7106 lines with **77** `g_bigstack`/`ThreadGroup`/CLI references and ~15 file-local static helpers. Pulling the solvers into a standalone `namespace plink2` TU that links only `plink2_matrix.cc` + pgenlib was the pragmatic way to reuse plink2's exact numerics without its thread-pool/CLI model.
- **Score:** 🟡 **YELLOW.** This is the "extracted-and-vendored" category by definition — correct and faithful, but a maintenance/drift burden. It should track the submodule pin (e.g. a comment recording the exact upstream line range + commit, and a re-sync check when the pin moves).

### 🟡 `vcf_genotype_parse.cpp` — vendored VCF GT parser

- **Feature:** biallelic VCF `GT`/dosage line parsing (unphased + phased) and quality-field checks.
- **Provenance:** header comment + inline banners state it copies `VcfConvertUnphasedBiallelicLine`, `VcfConvertPhasedBiallelicLine`, `VcfCheckQuals`, `VcfImportBaseContext`, and `VcfParseErr` from `plink2_import.cc`. Upstream functions are non-static (`plink2_import.cc:1480,1962`) but the file is a large CLI/import-pipeline TU; the copied functions "depend only on header-inline utilities from `plink2_base.h`/`plink2_string.h`."
- **Score:** 🟡 **YELLOW.** Same extracted-and-vendored drift burden as the GLM math TU; extraction is defensible because `plink2_import.cc` cannot be linked wholesale.

### 🟡 `plink_ld.cpp` — pairwise LD (r², D')

- **Feature:** r², D' and observation count for a variant pair or windowed pairs.
- **Wraps:** genotype fetch via `PgrGet` (`:551`); bit-layout constants (`DivUp`, `NypCtToAlignedWordCt`).
- **Reimplements:** `ComputeLdStats` (`:39–121`) walks 2-bit words inline, accumulates sums, and computes a **genotype-level Pearson** r² = cov²/(varₐ·var_b) (`:95`) plus a **Weir-1979 composite** D' (`:102–117`).
- **Callable upstream equivalent:** none usable. `plink2_ld.cc` entry points (`IndepPairwise`, `LdConsole`, …) take plink2-internal structs, use `g_bigstack`, `SetThreadFuncAndData`, and are pruning/console-oriented — not a per-pair r²/D' function. Not linked.
- **Score:** 🟡 **YELLOW.** Reimplementation is forced, but flagged for review: the genotype-covariance estimator differs from plink2's phased-haplotype r², and the code itself notes D' can exceed 1 under HWE deviation. Correctness-drift concern independent of wrappability.

### 🟡 `plink_pca.cpp` — PCA

- **Feature:** top-`n_pcs` eigenvectors/eigenvalues of the sample GRM via multi-pass randomized subspace iteration; per-sample PC output.
- **Wraps:** `PgrGetCounts` (`:375`, freq + monomorphic skip), `PgrGet` (`:878`), `GenoarrToBytesMinus9` (`:884`).
- **Reimplements:** genotype normalization (`ComputeVariantNorm`/`NormalizeGenotypes`, `plink_common.cpp:1432/1446`), a bespoke DuckDB-threaded randomized-subspace accumulator (`AccumulateStepA/B`, `MergePass`), and the SVD via **Eigen `BDCSVD`** (`RunKrylovSVD` `:667`, `RunFinalSVD` `:684`) — not LAPACK, not plink2.
- **Callable upstream equivalent:** none usable. `CalcPca` (`plink2_matrix_calc.cc`) is saturated with `g_bigstack_*`, `PgenMtLoadInit`, `SetThreadFuncAndData(...Thread)`; the file is not compiled. Approach *mirrors* plink2's `--pca approx` randomized SVD but the implementation is independent.
- **Score:** 🟡 **YELLOW.** Substantial independent numerical code (fails GREEN's "small + faithful") with genuine drift risk, but no callable equivalent exists → not RED. Recommend cross-checking eigenvalues/vectors against `plink2 --pca approx` on a fixture.

### 🟢 `pgen_reader.cpp` — genotype decode

- **Feature:** `read_pgen()` — decode `.pgen` records to hardcalls / dosages / phased pairs.
- **Wraps:** `PgrGet`, `PgrGetD` (`:844`, dosages), `PgrGetP` (phased), `GenoarrToBytesMinus9` (`:870/881`); full `Pgfi/Pgr` lifecycle; correct buffer sizing helpers. The only local transform is `UnpackPhasedGenotypes` (see plink_common) — a thin presentation reshape, not an algorithm.
- **Score:** 🟢 **GREEN.** Textbook thin wrapper over pgenlib decode.

### 🟢 `plink_freq.cpp` — allele frequency

- **Feature:** per-variant ALT frequency + counts.
- **Wraps:** `PgrGetCounts` / `PgrGetDCounts` (`:469/474`); AF is a trivial ratio over `genocounts` (`:534`), denom `2·non-missing` per plink semantics.
- **Score:** 🟢 **GREEN.** Counting is delegated; the arithmetic on the returned counts is unavoidable glue.

### 🟢 `plink_missing.cpp` — missingness

- **Feature:** per-variant and per-sample missingness / F_MISS.
- **Wraps:** `PgrGetMissingness` (`:` throughout) + `PopcountWords`; F_MISS is a division.
- **Score:** 🟢 **GREEN.** Delegates the hard part (bitarray missingness) to pgenlib.

### 🟢 `plink_score.cpp` — polygenic risk score

- **Feature:** per-sample Σ weightᵢ·dosageᵢ with flip / center / mean-imputation modes.
- **Wraps:** `PgrGetD` (`:573`), `Dosage16ToDoublesMinus9` (`:581`).
- **Reimplements:** the weighted dot-product accumulation + imputation branches (`:562–685`). `pgenlib_ffi_support.h` exposes `LinearCombinationMeanimpute`, but it computes a single per-variant combination, not the extension's per-sample-across-variants PRS with allele-flip/center columns; `plink2_score.cc` isn't even in the tree.
- **Score:** 🟢 **GREEN.** Genotype/dosage read is wrapped; the remaining code is a straightforward dot product — unavoidable glue with no callable per-sample-PRS equivalent.

### 🟢 `plink_glm.cpp` — GWAS regression (core, excluding the p-value layer noted above)

- **Wraps (genotype I/O):** `PgrGetD` (`:1310`), `Dosage16ToDoublesMinus9` (`:1317`), full `Pgfi/Pgr` lifecycle.
- **Wraps (linear algebra):** links `plink2_matrix.cc` and calls `MultiplySelfTranspose` (`:1092`), `LinearRegressionInv` (`:1101`), `InvertSymmdefFmatrixFirst/SecondHalf` (`:1247/1252`) for multivariate OLS + Hessian inverse.
- **Wraps (logistic/Firth):** calls `LogisticRegressionF` (`:1202`) then `FirthRegressionF` (`:1223`) from the vendored solver TU.
- **Minor reimplementation note:** the *no-covariate* univariate linear case is a hand-rolled closed form (`:1008–1050`) rather than routing through `LinearRegressionInv`. Trivial/textbook and low-risk; a small YELLOW-flavored optimization (could reuse the multivariate path) but not a correctness concern.
- **Score:** 🟢 **GREEN** for the regression engine itself — it wraps plink2's matrix lib and vendors plink2's exact solver. (The distribution layer is the separate 🔴 above.)

### 🟢 `pvar_reader.cpp` / `psam_reader.cpp` — text metadata readers

- **Feature:** `read_pvar()` / `read_psam()` (+ `.bim`/`.fam`) → DuckDB columns, including arbitrary INFO/FILTER/custom columns, typed per header, with a streaming path for ~10GB `.pvar`.
- **Reimplements:** tab/whitespace line splitting + header/type inference (`ParsePvarHeader`, `PvarColumnType`).
- **Callable upstream equivalent:** `pvar_ffi_support.h` exposes `LoadMinimalPvarEx`, but it loads only *minimal* variant info (ID/pos/alleles) into plink2 memory structures — it cannot surface arbitrary INFO/FILTER columns as DuckDB vectors. Upstream `plink2_pvar.cc`/`plink2_psam.cc` are CLI/bigstack-bound.
- **Score:** 🟢 **GREEN.** This is DuckDB-native CSV-style ingestion, not a duplication of a plink2 algorithm; the "reimplemented" logic is trivial tab-splitting driven by extension-specific schema needs.

### 🟢 `pfile_reader.cpp` — unified reader (3699 lines)

- **Feature:** `read_pfile()` — orient (variant/genotype/sample) + genotype-mode + filter orchestration over the pgen/pvar/psam trio.
- **Audit sweep:** grep for statistical markers (`sqrt|chisq|pval|erfc|lgamma|betai|normal|cdf|correl|hwe|fisher|regress`) returns **no** numeric-algorithm hits — it is pure I/O + reshape/filter glue, delegating all decode to the same pgenlib calls and `UnpackPhasedGenotypes`.
- **Score:** 🟢 **GREEN.** Largest file, but entirely marshaling/orchestration.

### 🟢 `plink_common.cpp` — shared P2 infrastructure (2256 lines)

- **Feature:** RAII buffers, sample-subset (`sample_include` + interleaved vec + cumulative popcounts), region filtering, companion-file discovery, `UnpackPhasedGenotypes` (`:1460`), sex-aware genotype-count classification (`:~2090`), and the PCA normalization helper `ComputeVariantNorm` (`:1441`, `1/sqrt(2·f·(1-f))`).
- **Assessment:** the only "math" is the trivial standardization constant and genotype-into-HWE-strata *counting* (not a statistic). `UnpackPhasedGenotypes` is a 35-line reshape of already-decoded bytes into the extension's `[allele1,allele2]` convention — distinct in I/O and output type from the linked `GenoarrPhasedToAlleleCodes` (which emits table-driven int32 codes), so adapting it would add work, not remove it.
- **Score:** 🟢 **GREEN.** Infrastructure + thin transforms built directly on pgenlib primitives.

### 🟢 `vcf_reader.cpp` — `read_plink_vcf()`

- **Feature:** stream a VCF, parse biallelic `GT` (and dosage) fields, and emit genotypes in the same array/list/columns modes as the pgen readers.
- **Wraps:** builds a `genovec` then decodes via `plink2::GenoarrToBytesMinus9` (`:693`) and reuses `UnpackPhasedGenotypes` (`:696`); includes `pgenlib_misc.h` + `pgenlib_ffi_support.h` and uses `NypCtToAlignedWordCt`/`BitCtToAlignedWordCt` for buffer sizing (`:433/437`). The actual GT field parsing is delegated to the vendored `vcf_genotype_parse.cpp` (scored 🟡 above).
- **Audit sweep:** the statistical-marker grep (`sqrt|chisq|pval|erfc|lgamma|betai|normal|cdf|correl|hwe|fisher|regress`) returns **no** hits — the only local logic is VCF text scanning (`ParseVcfRegion`, `ParseVcfPos`, tab-position indexing).
- **Score:** 🟢 **GREEN.** Wraps pgenlib for the genovec→bytes decode and reuses the vendored parser; no statistical reimplementation of its own. (Counts/stats modes are intentionally unsupported here precisely because they'd need `PgrGetCounts`, `:306`.)

### 🟢 `plinking_duck_extension.cpp` — entry point

- Registration of table functions + config options + an Eigen3 availability guard. No algorithm. 🟢 **GREEN.**

---

## Cross-cutting notes

1. **One build decision drives every RED.** Both RED findings (Hardy HWE, GLM p-values) exist
   solely because `plink2_stats.cc` is excluded (`CMakeLists.txt:128`). That file appears
   cleanly linkable (only `plink2_stats.h` + `plink2_string.h` + libc; 0 bigstack/thread refs).
   Linking it exposes callable `HweLnP`, `HweXchrLnP`, `FisherExact2x2P`, `TstatToP2`,
   `ChisqToP`, `QuantileToZscore` — enabling wraps that also *improve* correctness (overflow
   safety, proper chrX HWE). This is the single highest-value change in the audit.
   - **Caveat to weigh before linking:** `plink2_stats.cc` includes `<stdlib.h> // exit()` and
     some upstream paths call `exit()` on error rather than returning. A wrapper would need to
     validate inputs up front (or guard) so a bad call can't terminate the DuckDB process.
     Does not change any score, but is a real integration constraint.

2. **Vendored extracts are the accepted pattern, and it's a defensible one.**
   `plink2_glm_logistic_math.cpp` and `vcf_genotype_parse.cpp` both copy non-static upstream
   functions because their home TUs (`plink2_glm_logistic.cc` at 7106 lines / 77 bigstack refs;
   `plink2_import.cc` CLI pipeline) cannot be linked wholesale. This is the correct trade-off
   for plink2's architecture — but each is a manual-resync liability. Recommend: record the exact
   upstream file + line range + submodule commit at the top of each vendored TU, and add a
   pin-bump checklist item to diff against upstream.

3. **Forced reimplementations (PCA, LD) are fair YELLOWs, not REDs.** plink2's PCA
   (`CalcPca`) and LD (`plink2_ld.cc`) are inseparable from `g_bigstack` + `PgenMtLoad` +
   `SetThreadFuncAndData` + CLI. There is no callable equivalent, so hand-rolling is the only
   option short of vendoring those subsystems too. They are YELLOW (not GREEN) only because the
   reimplemented math is substantial and independent — carrying drift risk that warrants
   fixture-level cross-checks against plink2 output, and (for LD) a review of the
   genotype-covariance-vs-phased-haplotype r² choice.

4. **The genotype/dosage read path is uniformly GREEN.** Every reader and stat function decodes
   through pgenlib (`PgrGet*`, `GenoarrToBytesMinus9`, `Dosage16ToDoublesMinus9`) with correct
   buffer-sizing helpers. No file re-implements the 2-bit/dosage decode — the core wrapping
   discipline is intact.
