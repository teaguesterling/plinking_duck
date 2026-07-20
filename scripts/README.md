# scripts/

Maintenance scripts for the PlinkingDuck extension.

## `check_vendored_drift.sh` — vendored plink2 drift canary

Two source files are hand-copied ("vendored") extracts of upstream plink-ng (plink2)
functions, because their home translation units are huge and entangled with plink2's
`g_bigstack` arena / `ThreadGroup` / CLI model and cannot be linked wholesale (see
`docs/planning/plink2-wrapping-audit.md`, YELLOW #2/#3):

| Vendored file | Upstream source |
| --- | --- |
| `src/plink2_glm_logistic_math.cpp` | `third_party/plink-ng/2.0/plink2_glm_logistic.cc` |
| `src/vcf_genotype_parse.cpp` | `third_party/plink-ng/2.0/plink2_import.cc` |

The `third_party/plink-ng` submodule is pinned to a fixed commit, so upstream only changes
when someone bumps the submodule. This canary catches exactly that: a submodule bump that
alters any vendored function body.

Run it:

```sh
bash scripts/check_vendored_drift.sh
```

- Exit `0` with `OK: <file> matches pinned upstream` when every vendored function body is
  byte-identical to the pinned upstream (modulo documented, behavior-neutral renamings —
  see the script header and each vendored file's `VENDORED-CODE PIN` block).
- Exit `1` with a message naming the file and the drifted function(s), plus re-extract
  instructions, when they differ.

How it works: for each named function it locates the definition in both the vendored copy
and the current submodule source by a column-0 `<ret> Name(` marker (robust to line-number
shifts), captures the body from the first `{` to the top-level closing `}`, normalizes
(strips comments and all whitespace), and compares sha256s. See the script's header comment
for the full method, the two VCF renaming rules, and documented limitations (signatures and
header-level macro/inline-helper changes are out of scope).

When it fails after an intentional submodule bump, follow the **Re-sync procedure** in the
`VENDORED-CODE PIN` block at the top of the affected `.cpp`.

### CI wiring (optional)

This is a clean, dependency-light check (`bash` + `awk` + `perl` + `sha256sum`). To run it in
CI, add a step to a submodule-checked-out job in
`.github/workflows/MainDistributionPipeline.yml` before the build:

```yaml
- name: Check vendored plink2 drift
  run: bash scripts/check_vendored_drift.sh
```

It is documented here rather than auto-added so it can be slotted into the correct job
without perturbing the release matrix.
