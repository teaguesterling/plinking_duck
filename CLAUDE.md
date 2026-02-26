# PlinkingDuck

DuckDB extension for reading PLINK 2 genomics file formats in SQL.

## Project Structure

- `src/` — C++ extension source
  - `plinking_duck_extension.cpp` — entry point, registers all table functions
  - `pvar_reader.cpp` / `.hpp` — `read_pvar()`: .pvar and .bim files
  - `psam_reader.cpp` / `.hpp` — `read_psam()`: .psam and .fam files
  - `pgen_reader.cpp` / `.hpp` — `read_pgen()`: .pgen binary genotype files (uses pgenlib)
- `test/sql/` — DuckDB sqllogictest files (positive + negative per reader)
- `test/data/` — test fixtures (small VCF-derived PLINK files)
- `docs/planning/` — design docs and implementation plans
- `third_party/plink-ng/` — pgenlib submodule (read path only)

## Build & Test

- `make -j 4` to build (blq command: `build`)
- `make test` to run all sqllogictests (blq command: `test`)
- Build outputs: `./build/release/duckdb`, `./build/release/test/unittest`
- pgenlib_write.cc is intentionally excluded from the build

## Architecture Patterns

- **Table function lifecycle**: Bind → InitGlobal → InitLocal (per-thread) → Scan
- **Projection pushdown**: all readers support it; skip expensive work (e.g., PgrGet) when columns aren't referenced
- **VFS integration**: use DuckDB FileSystem for all file I/O (not raw stdio)
- **Bulk file reads**: read entire file into memory then parse (see ReadFileLines in psam_reader, pgen_reader)
- **Parallel scan**: atomic batch claiming on variant index; per-thread local state

## pgenlib Notes

- PgenFileInfo must outlive PgenReader (pgr holds reference to pgfi)
- Per-thread: each thread needs its own PgenFileInfo + PgenReader (FILE* not thread-safe)
- Buffer sizing: use `NypCtToAlignedWordCt` (not `DivUp`) for SIMD-safe genovec allocation
- Cleanup order: CleanupPgr before CleanupPgfi
- `GenoarrToBytesMinus9`: unpacks 2-bit genovec → int8 array (-9 = missing)

## Test Data

- `pgen_example.*`: 4 variants × 4 samples (main pgen test dataset)
- `pgen_orphan.*`: .pgen + .pvar only, no .psam (index-only mode)
- `all_missing.*`: 2 variants × 2 samples, all genotypes missing
- `example.*`: pvar/psam/bim/fam test data for P1-001/P1-002 — do not overwrite
- `generate_test_data.sh`: regenerates pgen fixtures (needs plink2)
- `plink2`: Binary is in `/mnt/aux-data/teague/Dev/spack/opt/spack/linux-zen3/plink2-2.0.0-a.6.9-e2l3wx22vy6tkmdwzhaexlml4rs3okx3`

## Conventions

- sqllogictest format specifiers: `T` = varchar, `I` = integer, `R` = real
- DuckDB function overloads: can't register two with same positional args; use `LogicalType::ANY` for type-dispatched named params
- Named parameters go in `read_pgen.named_parameters["key"]`

<!-- blq:agent-instructions -->
## blq - Build Log Query

Run builds and tests via blq MCP tools, not via Bash directly:
- `mcp__blq_mcp__commands` - list available commands
- `mcp__blq_mcp__run` - run a registered command (e.g., `run(command="test")`)
- `mcp__blq_mcp__register_command` - register new commands
- `mcp__blq_mcp__status` - check current build/test status
- `mcp__blq_mcp__errors` - view errors from runs
- `mcp__blq_mcp__info` - detailed run info (supports relative refs like `+1`, `latest`)
- `mcp__blq_mcp__output` - search/filter captured logs (grep, tail, head, lines)

Do NOT use shell pipes or redirects in commands (e.g., `pytest | tail -20`).
Instead: run the command, then use `output(run_id=N, tail=20)` to filter.
<!-- /blq:agent-instructions -->
