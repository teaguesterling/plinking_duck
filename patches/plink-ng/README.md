# plink-ng vendored patches

Minimal, re-appliable patches to the pinned `third_party/plink-ng` submodule. We
keep pgenlib pristine except for these; each is meant to be **upstreamable**.

## Mechanism

`CMakeLists.txt` applies every `patches/plink-ng/*.patch` (sorted) at **configure**
time, idempotently: `git apply --reverse --check` succeeds only when a patch is
already present, so it is applied exactly once — safe on reconfigure, a fresh
`git submodule update`, and CI from a clean checkout. The submodule pointer is
**not** modified in the parent repo; the patch is the durable record and the build
applies it to the working tree.

To iterate on a patch: edit the submodule file in place, rebuild/test, then
regenerate the patch with
`git -C third_party/plink-ng diff -- <file> > patches/plink-ng/NNNN-name.patch`
and `git -C third_party/plink-ng checkout -- <file>` (CMake re-applies on next
configure).

## Patches

- **0001-pgen-vfs-fopen-hook.patch** — routes pgenlib_read.cc's 3 `.pgen`/`.pgi`
  `fopen()` sites through `plinking_pgen_fopen()` (a declaration + 3 one-line
  redirects). Default behavior is unchanged (real `fopen`); when the extension
  registers a per-thread VFS opener, opens are served by a `fopencookie`/`funopen`
  `FILE*` over a DuckDB `FileHandle` (remote/VFS `.pgen` reads). The hook is defined
  in `src/plinking_pgen_vfs.{hpp,cpp}` (compiled into the `pgenlib` static lib, no
  DuckDB dependency). See `docs/planning/P12-001-*`.

## Drift canary

`scripts/check_vendored_drift.sh` verifies our vendored *extracts* against the
submodule; these patches are the intended delta to the submodule itself. A submodule
bump that moves the patched `fopen` lines will make `git apply` fail loudly at
configure — regenerate the patch against the new revision.
