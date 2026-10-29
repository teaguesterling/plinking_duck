# cmake-format configuration for plinking_duck.
#
# A .py config (not .yaml): cmakelang loads a YAML config by importing PyYAML, and
# DuckDB's format venv does not install it -- duckdb/Makefile builds format_venv with
# black, cmake-format and clang_format only, and format.py treats ANY formatter stderr
# as fatal. A .cmake-format.yaml therefore CRASHES the Format Check job rather than
# failing it. JSON and Python configs need no PyYAML.
#
# enable_markup is OFF deliberately. With it on, cmake-format treats comment bodies as
# markup and REFLOWS them, which ran the numbered list of the three plink-ng static
# libraries together into one paragraph and broke sentences mid-clause. The comments in
# CMakeLists.txt carry measured reasons, not decoration.
#
# tab_size is 4 to match the indentation this file already uses. Measured: with the
# default of 2 the formatter rewrote 270 lines, 20 of them comment lines; with 4 it
# rewrites 169 and leaves every comment byte-for-byte identical.
#
# Layout of the CMake code itself is still checked and enforced.
enable_markup = False
line_width = 100
tab_size = 4
