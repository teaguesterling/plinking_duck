PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Configuration of extension
EXT_NAME=plinking_duck
EXT_CONFIG=${PROJ_DIR}extension_config.cmake

# Include the Makefile from extension-ci-tools
include extension-ci-tools/makefiles/duckdb_extension.Makefile

# ---------------------------------------------------------------------------
# tidy-check override (real coverage).
#
# The upstream tidy-check (extension-ci-tools) selects translation units with the
# regex 'src/.*/'. Its trailing slash requires a SUBDIRECTORY of src/, but
# plinking's sources are flat (src/*.cpp), so that regex matches zero TUs and the
# gate lints nothing while still reporting green. This override:
#   - selects the flat sources with 'src/[^/]+\.cpp$';
#   - anchors the header filter to plinking's own src (excludes duckdb-core headers);
#   - forces -std=c++17 (the tidy build configures C++14, so without this clang-tidy
#     emits spurious clang-diagnostic-c++17-extensions that would fail the gate);
#   - copies the repo-root .clang-tidy (now a real, owned file, not the duckdb symlink);
#   - PRINTS the TUs it selected and FAILS if the count is zero, so a green run can
#     never again be a silently empty gate.
tidy-check:
	mkdir -p ./build/tidy
	cmake $(GENERATOR) $(BUILD_FLAGS) $(EXT_DEBUG_FLAGS) -DDISABLE_UNITY=1 -DCLANG_TIDY=1 -DCMAKE_CXX_STANDARD=17 -S $(DUCKDB_SRCDIR) -B build/tidy
	cp .clang-tidy build/tidy/.clang-tidy
	@python3 -c "import json,re,sys; db=json.load(open('build/tidy/compile_commands.json')); rx=re.compile(r'$(PROJ_DIR)src/[^/]+[.]cpp$$'); fs=sorted({e['file'] for e in db if rx.search(e['file'])}); print('[tidy] translation units selected:', len(fs)); [print('   ', f) for f in fs]; sys.exit('ERROR: tidy selected 0 translation units — the source regex matched nothing (hollow gate)' if not fs else 0)"
	cd build/tidy && python3 ../../duckdb/scripts/run-clang-tidy.py '$(PROJ_DIR)src/[^/]+\.cpp$$' -header-filter '$(PROJ_DIR)src/.*' -extra-arg=-std=c++17 -quiet $(TIDY_THREAD_PARAMETER) $(TIDY_BINARY_PARAMETER) $(TIDY_PERFORM_CHECKS)