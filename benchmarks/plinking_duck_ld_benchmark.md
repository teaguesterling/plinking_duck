## Purpose

This benchmark compares local builds of:

- upstream `teaguesterling/plinking_duck` cloned at
  `.sync/plinking_duck`
- `sounkou-bioinfo/plinking_duck` fork cloned at
  `.sync/plinking_duck_sounkou`

The fork adds signed `R`, projection-aware LD output, and a bounded
cached popcount path for no-missing regional LD scans. All timings below
are measured with `/usr/bin/time -v`, so the table includes wall time
and peak resident set size.

## Machine, inputs, and build artifacts

``` r
threads <- as.integer(Sys.getenv("PLINKING_BENCH_THREADS", "19"))
n_reps <- as.integer(Sys.getenv("PLINKING_BENCH_REPS", "3"))

paths <- list(
  upstream_dir = normalizePath(".sync/plinking_duck"),
  fork_dir = normalizePath(".sync/plinking_duck_sounkou"),
  upstream_duckdb = normalizePath(".sync/plinking_duck/build/release/duckdb"),
  fork_duckdb = normalizePath(".sync/plinking_duck_sounkou/build/release/duckdb"),
  pgen = normalizePath(".sync/plink_resources/1kg_chr22_hg38/chr22_hg38.pgen"),
  pvar = normalizePath(".sync/plink_resources/1kg_chr22_hg38/chr22_hg38.pvar"),
  psam = normalizePath(".sync/plink_resources/1kg_chr22_hg38/chr22_hg38.psam"),
  plink2 = normalizePath(".sync/plink2/plink2", mustWork = FALSE)
)

stopifnot(file.exists(paths$upstream_duckdb), file.exists(paths$fork_duckdb),
          file.exists(paths$pgen), file.exists(paths$pvar), file.exists(paths$psam))

cmd_out <- function(cmd, args = character()) {
  paste(system2(cmd, args, stdout = TRUE, stderr = TRUE), collapse = "\n")
}

git_head <- function(dir) cmd_out("git", c("-C", dir, "rev-parse", "--short", "HEAD"))
file_mb <- function(path) round(file.info(path)$size / 1024^2, 2)

info <- data.frame(
  item = c("detected CPUs", "benchmark threads", "repetitions", "upstream HEAD", "fork HEAD",
           "upstream duckdb MB", "fork duckdb MB", "PGEN MB", "region", "window kb", "r2 threshold"),
  value = c(cmd_out("nproc"), threads, n_reps, git_head(paths$upstream_dir), git_head(paths$fork_dir),
            file_mb(paths$upstream_duckdb), file_mb(paths$fork_duckdb), file_mb(paths$pgen),
            "22:11300000-11400000", 100, 0.2)
)
knitr::kable(info)
```

| item               | value                |
|:-------------------|:---------------------|
| detected CPUs      | 20                   |
| benchmark threads  | 19                   |
| repetitions        | 3                    |
| upstream HEAD      | 49eab7b              |
| fork HEAD          | cd95239              |
| upstream duckdb MB | 48.4                 |
| fork duckdb MB     | 48.45                |
| PGEN MB            | 126.42               |
| region             | 22:11300000-11400000 |
| window kb          | 100                  |
| r2 threshold       | 0.2                  |

## pgenlibr header/linkage audit

This answers the linkage question directly: the installed R package
exposes R functions and a shared object, but it does **not** install
public C/C++ headers under `include/`. The CRAN source tarball contains
PLINK2/pgenlib headers under `tools/include`, but those are source-build
internals rather than a stable `LinkingTo: pgenlibr` API. For a DuckDB
extension, linking directly to vendored PLINK2 `pgenlib` is the right
boundary.

``` r
pgenlibr_installed <- requireNamespace("pgenlibr", quietly = TRUE)
pgenlibr_pkg_path <- if (pgenlibr_installed) system.file(package = "pgenlibr") else NA_character_
pgenlibr_installed_include <- if (pgenlibr_installed) file.path(pgenlibr_pkg_path, "include") else NA_character_
pgenlibr_source_headers <- list.files(".sync/pgenlibr", pattern = "\\.(h|hpp)$", recursive = TRUE, full.names = FALSE)
pgenlib_source_headers <- grep("pgenlib|plink2|pvar", pgenlibr_source_headers, value = TRUE)

audit <- data.frame(
  check = c("pgenlibr installed", "installed package path", "installed include/ exists",
            "source headers matching pgenlib/plink2/pvar"),
  value = c(pgenlibr_installed, pgenlibr_pkg_path, dir.exists(pgenlibr_installed_include),
            length(pgenlib_source_headers))
)
knitr::kable(audit)
```

| check                                       | value                            |
|:--------------------------------------------|:---------------------------------|
| pgenlibr installed                          | TRUE                             |
| installed package path                      | /usr/lib/R/site-library/pgenlibr |
| installed include/ exists                   | FALSE                            |
| source headers matching pgenlib/plink2/pvar | 15                               |

``` r
head(pgenlib_source_headers, 20)
```

    ##  [1] "src/pvar.h"                                  "tools/include/include/pgenlib_ffi_support.h"
    ##  [3] "tools/include/include/pgenlib_misc.h"        "tools/include/include/pgenlib_read.h"       
    ##  [5] "tools/include/include/plink2_base.h"         "tools/include/include/plink2_bgzf.h"        
    ##  [7] "tools/include/include/plink2_bits.h"         "tools/include/include/plink2_float.h"       
    ##  [9] "tools/include/include/plink2_htable.h"       "tools/include/include/plink2_memory.h"      
    ## [11] "tools/include/include/plink2_string.h"       "tools/include/include/plink2_text.h"        
    ## [13] "tools/include/include/plink2_thread.h"       "tools/include/include/plink2_zstfile.h"     
    ## [15] "tools/include/include/pvar_ffi_support.h"

## Benchmark helpers

``` r
parse_elapsed <- function(x) {
  if (length(x) == 0 || is.na(x)) return(NA_real_)
  parts <- as.numeric(strsplit(x, ":", fixed = TRUE)[[1]])
  if (length(parts) == 3) return(parts[1] * 3600 + parts[2] * 60 + parts[3])
  if (length(parts) == 2) return(parts[1] * 60 + parts[2])
  parts[1]
}

field <- function(lines, prefix) {
  hit <- grep(paste0("^\\s*", prefix), lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(paste0("^\\s*", prefix, ":\\s*"), "", hit[1])
}

run_timed <- function(label, exe, args, rep = NA_integer_) {
  stdout_file <- tempfile("bench_stdout_")
  stderr_file <- tempfile("bench_stderr_")
  # system2() still invokes a shell on Unix, so quote every executable/argument;
  # SQL strings contain semicolons and parentheses.
  status <- system2("/usr/bin/time", c("-v", shQuote(exe), shQuote(args)), stdout = stdout_file, stderr = stderr_file)
  if (is.null(status)) status <- 0L
  out <- if (file.exists(stdout_file)) readLines(stdout_file, warn = FALSE) else character()
  err <- if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character()

  elapsed_raw <- field(err, "Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\)")
  data.frame(
    label = label,
    rep = rep,
    status = as.integer(status),
    elapsed_sec = parse_elapsed(elapsed_raw),
    user_sec = as.numeric(field(err, "User time \\(seconds\\)")),
    sys_sec = as.numeric(field(err, "System time \\(seconds\\)")),
    cpu_percent = as.numeric(sub("%", "", field(err, "Percent of CPU this job got"), fixed = TRUE)),
    max_rss_mb = as.numeric(field(err, "Maximum resident set size \\(kbytes\\)")) / 1024,
    stdout = paste(out, collapse = "\\n"),
    stderr_tail = paste(unique(c(grep("Error|Exception|Binder", err, value = TRUE), tail(err, 12))), collapse = "\\n"),
    stringsAsFactors = FALSE
  )
}

ld_from <- function(extra = "") {
  sprintf(
    "FROM plink_ld('%s', pvar := '%s', psam := '%s', region := '22:11300000-11400000', window_kb := 100, r2_threshold := 0.2%s)",
    paths$pgen, paths$pvar, paths$psam, extra
  )
}

duck_sql <- function(select, extra = "") {
  paste0("SET threads=", threads, "; SET plinking_max_threads=", threads, "; ", select, " ", ld_from(extra), ";")
}

run_duckdb <- function(repo, exe, scenario, select, extra = "", reps = n_reps) {
  rows <- vector("list", reps)
  for (i in seq_len(reps)) {
    sql <- duck_sql(select, extra)
    rows[[i]] <- run_timed(paste(repo, scenario, sep = ":"), exe, c("-csv", "-c", sql), rep = i)
    rows[[i]]$repo <- repo
    rows[[i]]$scenario <- scenario
  }
  do.call(rbind, rows)
}

summarise_bench <- function(x) {
  ok <- x[x$status == 0, ]
  parts <- split(ok, paste(ok$scenario, ok$repo, sep = "\r"))
  out <- lapply(parts, function(d) {
    data.frame(
      scenario = d$scenario[1],
      repo = d$repo[1],
      n = nrow(d),
      median_elapsed_sec = median(d$elapsed_sec),
      min_elapsed_sec = min(d$elapsed_sec),
      median_max_rss_mb = median(d$max_rss_mb),
      median_cpu_percent = median(d$cpu_percent),
      output = gsub("\\n", "; ", d$stdout[1], fixed = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
```

## Run DuckDB LD benchmarks

Scenarios:

- `count_pairs`: compute all windowed LD pairs and return only the
  count.
- `sum_r2`: compute and aggregate `R2` only.
- `sum_r2_dprime`: compute `R2` plus `D_PRIME`, validating that the
  optimized path preserves the legacy statistic when requested.

``` r
benchmarks <- list(
  run_duckdb("upstream", paths$upstream_duckdb, "count_pairs", "SELECT COUNT(*)"),
  run_duckdb("fork", paths$fork_duckdb, "count_pairs", "SELECT COUNT(*)"),
  run_duckdb("upstream", paths$upstream_duckdb, "sum_r2", "SELECT SUM(R2)"),
  run_duckdb("fork", paths$fork_duckdb, "sum_r2", "SELECT SUM(R2)"),
  run_duckdb("upstream", paths$upstream_duckdb, "sum_r2_dprime", "SELECT COUNT(*), SUM(R2), SUM(D_PRIME)"),
  run_duckdb("fork", paths$fork_duckdb, "sum_r2_dprime", "SELECT COUNT(*), SUM(R2), SUM(D_PRIME)")
)
raw <- do.call(rbind, benchmarks)
raw_out <- file.path("benchmarks", "plinking_duck_ld_benchmark_raw.csv")
write.csv(raw, raw_out, row.names = FALSE)
knitr::kable(raw[, c("scenario", "repo", "rep", "status", "elapsed_sec", "user_sec", "sys_sec", "cpu_percent", "max_rss_mb")], digits = 3)
```

| scenario      | repo     | rep | status | elapsed_sec | user_sec | sys_sec | cpu_percent | max_rss_mb |
|:--------------|:---------|----:|-------:|------------:|---------:|--------:|------------:|-----------:|
| count_pairs   | upstream |   1 |      0 |        5.59 |   102.45 |    0.36 |        1838 |    348.078 |
| count_pairs   | upstream |   2 |      0 |        6.85 |   125.06 |    0.35 |        1829 |    348.891 |
| count_pairs   | upstream |   3 |      0 |        6.76 |   124.85 |    0.35 |        1849 |    348.332 |
| count_pairs   | fork     |   1 |      0 |        1.55 |    28.95 |    0.04 |        1860 |    120.562 |
| count_pairs   | fork     |   2 |      0 |        1.55 |    28.95 |    0.04 |        1863 |    120.551 |
| count_pairs   | fork     |   3 |      0 |        1.55 |    28.86 |    0.05 |        1854 |    120.531 |
| sum_r2        | upstream |   1 |      0 |        6.77 |   125.02 |    0.34 |        1849 |    348.094 |
| sum_r2        | upstream |   2 |      0 |        6.89 |   127.08 |    0.34 |        1847 |    348.867 |
| sum_r2        | upstream |   3 |      0 |        6.86 |   126.41 |    0.35 |        1847 |    348.391 |
| sum_r2        | fork     |   1 |      0 |        1.57 |    29.27 |    0.04 |        1856 |    120.816 |
| sum_r2        | fork     |   2 |      0 |        1.58 |    29.47 |    0.05 |        1860 |    120.707 |
| sum_r2        | fork     |   3 |      0 |        1.58 |    29.29 |    0.06 |        1855 |    120.613 |
| sum_r2_dprime | upstream |   1 |      0 |        6.87 |   126.65 |    0.38 |        1846 |    348.176 |
| sum_r2_dprime | upstream |   2 |      0 |        6.94 |   127.72 |    0.39 |        1843 |    348.477 |
| sum_r2_dprime | upstream |   3 |      0 |        7.04 |   129.19 |    0.41 |        1838 |    347.746 |
| sum_r2_dprime | fork     |   1 |      0 |        1.59 |    29.34 |    0.05 |        1843 |    121.031 |
| sum_r2_dprime | fork     |   2 |      0 |        1.59 |    29.55 |    0.04 |        1859 |    121.145 |
| sum_r2_dprime | fork     |   3 |      0 |        1.56 |    29.02 |    0.06 |        1853 |    120.977 |

## Summary and speedup

``` r
summary <- summarise_bench(raw)
summary_out <- file.path("benchmarks", "plinking_duck_ld_benchmark_summary.csv")
write.csv(summary, summary_out, row.names = FALSE)
knitr::kable(summary[, c("scenario", "repo", "n", "median_elapsed_sec", "min_elapsed_sec", "median_max_rss_mb", "median_cpu_percent")], digits = 3)
```

|                       | scenario      | repo     |   n | median_elapsed_sec | min_elapsed_sec | median_max_rss_mb | median_cpu_percent |
|:----------------------|:--------------|:---------|----:|-------------------:|----------------:|------------------:|-------------------:|
| count_pairsfork       | count_pairs   | fork     |   3 |               1.55 |            1.55 |           120.551 |               1860 |
| count_pairsupstream   | count_pairs   | upstream |   3 |               6.76 |            5.59 |           348.332 |               1838 |
| sum_r2fork            | sum_r2        | fork     |   3 |               1.58 |            1.57 |           120.707 |               1856 |
| sum_r2upstream        | sum_r2        | upstream |   3 |               6.86 |            6.77 |           348.391 |               1847 |
| sum_r2_dprimefork     | sum_r2_dprime | fork     |   3 |               1.59 |            1.56 |           121.031 |               1853 |
| sum_r2_dprimeupstream | sum_r2_dprime | upstream |   3 |               6.94 |            6.87 |           348.176 |               1843 |

``` r
up <- summary[summary$repo == "upstream", c("scenario", "median_elapsed_sec", "median_max_rss_mb")]
fk <- summary[summary$repo == "fork", c("scenario", "median_elapsed_sec", "median_max_rss_mb")]
names(up) <- c("scenario", "upstream_elapsed", "upstream_rss")
names(fk) <- c("scenario", "fork_elapsed", "fork_rss")
compare <- merge(up, fk, by = "scenario")
compare$speedup_x <- compare$upstream_elapsed / compare$fork_elapsed
compare$rss_delta_mb <- compare$fork_rss - compare$upstream_rss
knitr::kable(compare, digits = 3)
```

| scenario      | upstream_elapsed | upstream_rss | fork_elapsed | fork_rss | speedup_x | rss_delta_mb |
|:--------------|-----------------:|-------------:|-------------:|---------:|----------:|-------------:|
| count_pairs   |             6.76 |      348.332 |         1.55 |  120.551 |     4.361 |     -227.781 |
| sum_r2        |             6.86 |      348.391 |         1.58 |  120.707 |     4.342 |     -227.684 |
| sum_r2_dprime |             6.94 |      348.176 |         1.59 |  121.031 |     4.365 |     -227.145 |

## Full all_hg38 1000G benchmark (real large resource)

This run uses the real full-resource files at repository root:

- `./all_hg38.pgen.zst` as the source PGEN (decompressed to
  `./all_hg38.pgen` because pgenlib rejects `.pgen.zst` magic bytes)
- `./all_hg38_rs.pvar.zst` as the source PVAR, with
  `./all_hg38_rs.slim.pvar` derived by retaining the real
  `#CHROM POS ID REF ALT` columns and dropping the huge INFO field
- `./hg38_corrected.psam` as sample metadata

The benchmark interval is `22:11300000-16300000` (~5 Mb; 91,882
variants), with a 100 kb LD window and `r2_threshold = 0.2`. This is
intentionally much larger than the earlier 100 kb smoke benchmark and
exercises the full all-hg38 PGEN/PVAR metadata scale. The CSV below is
produced from the `/usr/bin/time -v` run; re-running it takes several
minutes and requires the decompressed 8.9 GB PGEN plus the 2.0 GB slim
PVAR.

Commands used to prepare real inputs (not committed because they are
large):

``` sh
zstd -T0 -d -k all_hg38.pgen.zst -o all_hg38.pgen
zstd -dc -T0 all_hg38_rs.pvar.zst | \
  awk 'BEGIN{FS=OFS="\t"} /^##/ {next} /^#CHROM/ {print "#CHROM","POS","ID","REF","ALT"; next} !/^#/ {print $1,$2,$3,$4,$5}' \
  > all_hg38_rs.slim.pvar
```

``` r
all_hg38 <- read.csv(file.path("benchmarks", "plinking_duck_ld_benchmark_all_hg38.csv"))
knitr::kable(all_hg38, digits = 3)
```

| repo     | local_head | scenario           | source_pgen_zst   | pgen          | pvar                  | psam                | region               | variant_ct | window_kb | r2_threshold | n_pairs | sum_r2 | elapsed_sec | user_sec | sys_sec | cpu_percent | max_rss_mb |
|:---------|:-----------|:-------------------|:------------------|:--------------|:----------------------|:--------------------|:---------------------|-----------:|----------:|-------------:|--------:|-------:|------------:|---------:|--------:|------------:|-----------:|
| upstream | 49eab7b    | all_hg38_chr22_5mb | all_hg38.pgen.zst | all_hg38.pgen | all_hg38_rs.slim.pvar | hg38_corrected.psam | 22:11300000-16300000 |      91882 |       100 |          0.2 | 1270300 | 709755 |      149.93 |  2552.85 |   23.13 |        1718 |  21867.734 |
| fork     | cd95239    | all_hg38_chr22_5mb | all_hg38.pgen.zst | all_hg38.pgen | all_hg38_rs.slim.pvar | hg38_corrected.psam | 22:11300000-16300000 |      91882 |       100 |          0.2 | 1270300 | 709755 |      138.93 |  2533.18 |    4.65 |        1826 |    699.164 |

``` r
all_hg38_compare <- merge(
  all_hg38[all_hg38$repo == "upstream", c("scenario", "elapsed_sec", "max_rss_mb")],
  all_hg38[all_hg38$repo == "fork", c("scenario", "elapsed_sec", "max_rss_mb")],
  by = "scenario", suffixes = c("_upstream", "_fork")
)
all_hg38_compare$speedup_x <- all_hg38_compare$elapsed_sec_upstream / all_hg38_compare$elapsed_sec_fork
all_hg38_compare$rss_delta_mb <- all_hg38_compare$max_rss_mb_fork - all_hg38_compare$max_rss_mb_upstream
knitr::kable(all_hg38_compare, digits = 3)
```

| scenario           | elapsed_sec_upstream | max_rss_mb_upstream | elapsed_sec_fork | max_rss_mb_fork | speedup_x | rss_delta_mb |
|:-------------------|---------------------:|--------------------:|-----------------:|----------------:|----------:|-------------:|
| all_hg38_chr22_5mb |               149.93 |            21867.73 |           138.93 |         699.164 |     1.079 |    -21168.57 |

Note: the fork now uses native PVAR region pushdown for `plink_ld()`
when a `region` is supplied and shares immutable pgenlib metadata across
workers. On the real all-hg38 resource this cuts peak RSS substantially;
elapsed time is now dominated by the LD kernel / PGEN reads rather than
all-PVAR or per-worker pgen metadata materialization.

## Output agreement

``` r
outputs <- summary[, c("scenario", "repo", "output")]
knitr::kable(outputs)
```

|                       | scenario      | repo     | output                                                                         |
|:----------------------|:--------------|:---------|:-------------------------------------------------------------------------------|
| count_pairsfork       | count_pairs   | fork     | count_star(); 42258                                                            |
| count_pairsupstream   | count_pairs   | upstream | count_star(); 42258                                                            |
| sum_r2fork            | sum_r2        | fork     | sum(R2); 26135.51331225652                                                     |
| sum_r2upstream        | sum_r2        | upstream | sum(R2); 26135.513312256506                                                    |
| sum_r2_dprimefork     | sum_r2_dprime | fork     | count_star(),sum(R2),sum(D_PRIME); 42258,26135.51331225653,19214.231785830212  |
| sum_r2_dprimeupstream | sum_r2_dprime | upstream | count_star(),sum(R2),sum(D_PRIME); 42258,26135.513312256517,19214.231785830172 |

## Signed R smoke test

The fork should expose signed `R`; upstream should fail because the
column does not exist.

``` r
signed_sql <- duck_sql("SELECT R, R2", "")
signed_sql <- sub(";$", " LIMIT 5;", signed_sql)
signed_upstream <- run_timed("upstream:signed_R", paths$upstream_duckdb, c("-csv", "-c", signed_sql), rep = 1)
signed_fork <- run_timed("fork:signed_R", paths$fork_duckdb, c("-csv", "-c", signed_sql), rep = 1)
signed <- rbind(
  transform(signed_upstream, repo = "upstream"),
  transform(signed_fork, repo = "fork")
)
knitr::kable(signed[, c("repo", "status", "elapsed_sec", "max_rss_mb", "stdout", "stderr_tail")], digits = 3)
```

| repo     | status | elapsed_sec | max_rss_mb | stdout                                                                                              | stderr_tail                                                                                                                                                                                                                                                                                                                                                               |
|:---------|-------:|------------:|-----------:|:----------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| upstream |      1 |        0.16 |    204.992 |                                                                                                     | Binder Error: Referenced column “R” not found in FROM clause!Major (requiring I/O) page faults: 0Minor (reclaiming a frame) page faults: 48793Voluntary context switches: 72Involuntary context switches: 8Swaps: 0File system inputs: 0File system outputs: 8Socket messages sent: 0Socket messages received: 0Signals delivered: 0Page size (bytes): 4096Exit status: 1 |
| fork     |      0 |        1.59 |     47.027 | R,R2,0.6664583333333333,0.6664583333333333,0.6664583333333333,0.6664583333333333,0.6664583333333333 | Major (requiring I/O) page faults: 0Minor (reclaiming a frame) page faults: 9914Voluntary context switches: 116Involuntary context switches: 6Swaps: 0File system inputs: 0File system outputs: 8Socket messages sent: 0Socket messages received: 0Signals delivered: 0Page size (bytes): 4096Exit status: 0                                                              |

## Population-weighted LD smoke benchmark

This exercises the GAUSS-like population-weighted path added to the
fork. The `.psam` file has a `SuperPop` column; each superpopulation
below contributes 20% of the total sample weight regardless of how many
samples it has in the panel.

``` r
weighted_sql <- duck_sql(
  "SELECT COUNT(*) AS n_pairs, SUM(R2) AS sum_r2, SUM(R) AS sum_r",
  ", population_column := 'SuperPop', population_weights := 'AFR=0.20,AMR=0.20,EAS=0.20,EUR=0.20,SAS=0.20'"
)
weighted_fork <- run_timed("fork:population_weighted", paths$fork_duckdb, c("-csv", "-c", weighted_sql), rep = 1)
weighted_display <- weighted_fork
weighted_display$stdout <- gsub("\\n", "; ", weighted_display$stdout, fixed = TRUE)
knitr::kable(weighted_display[, c("label", "status", "elapsed_sec", "user_sec", "sys_sec", "cpu_percent", "max_rss_mb", "stdout")], digits = 3)
```

| label                    | status | elapsed_sec | user_sec | sys_sec | cpu_percent | max_rss_mb | stdout                                                            |
|:-------------------------|-------:|------------:|---------:|--------:|------------:|-----------:|:------------------------------------------------------------------|
| fork:population_weighted |      0 |        9.02 |   170.07 |    0.08 |        1884 |    120.945 | n_pairs,sum_r2,sum_r; 41353,26317.251282885896,32140.667690079747 |

## PLINK2 CLI reference point

This is not a DuckDB implementation, but it is a useful external
reference for how fast PLINK2’s native LD kernel is on the same input
and region.

``` r
run_plink2 <- function(reps = n_reps) {
  if (!file.exists(paths$plink2)) return(data.frame())
  prefix <- sub("\\.pgen$", "", paths$pgen)
  rows <- vector("list", reps)
  for (i in seq_len(reps)) {
    out_prefix <- tempfile("plink2_ld_")
    args <- c(
      "--pfile", prefix,
      "--chr", "22", "--from-bp", "11300000", "--to-bp", "11400000",
      "--ld-window-kb", "100", "--ld-window-r2", "0.2",
      "--threads", as.character(threads),
      "--r-unphased", "ref-based", "cols=id,ref,alt",
      "--out", out_prefix
    )
    rows[[i]] <- run_timed("plink2:ld", paths$plink2, args, rep = i)
    rows[[i]]$repo <- "plink2"
    rows[[i]]$scenario <- "plink2_r_unphased"
  }
  do.call(rbind, rows)
}
plink2_raw <- run_plink2()
if (nrow(plink2_raw)) {
  write.csv(plink2_raw, file.path("benchmarks", "plinking_duck_ld_benchmark_plink2.csv"), row.names = FALSE)
  knitr::kable(plink2_raw[, c("scenario", "repo", "rep", "status", "elapsed_sec", "user_sec", "sys_sec", "cpu_percent", "max_rss_mb")], digits = 3)
} else {
  cat("PLINK2 binary not found; skipped.\n")
}
```

| scenario          | repo   | rep | status | elapsed_sec | user_sec | sys_sec | cpu_percent | max_rss_mb |
|:------------------|:-------|----:|-------:|------------:|---------:|--------:|------------:|-----------:|
| plink2_r_unphased | plink2 |   1 |      0 |        0.17 |     0.66 |    0.07 |         419 |    238.594 |
| plink2_r_unphased | plink2 |   2 |      0 |        0.17 |     0.65 |    0.07 |         410 |    237.031 |
| plink2_r_unphased | plink2 |   3 |      0 |        0.18 |     0.66 |    0.08 |         402 |    239.062 |

## Notes

- The fork trades memory for speed on bounded regional scans by caching
  packed genotypes per worker. The benchmark’s RSS delta is therefore
  expected.
- The cached fast path is only used for variants with no missing
  genotypes; the code falls back to the legacy pairwise scanner when
  missingness is present.
- Dense LD for ColocBoost/SuSiE should still use block/tiled matrix
  kernels in a future pass; this benchmark covers the current
  `plink_ld()` pair-table API.
