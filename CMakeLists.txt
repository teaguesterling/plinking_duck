cmake_minimum_required(VERSION 3.5)

# Set extension name here
set(TARGET_NAME plinking_duck)

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})
include_directories(src/include)

set(EXTENSION_SOURCES
    src/plinking_duck_extension.cpp
    src/pvar_reader.cpp
    src/psam_reader.cpp
    src/pgen_reader.cpp
    src/pfile_reader.cpp
    src/plink_common.cpp
    src/plink_freq.cpp
    src/plink_hardy.cpp
    src/plink_missing.cpp
    src/plink_ld.cpp
    src/plink_score.cpp
)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

# ===========================================================================
# pgenlib build integration
#
# Three static libraries are built from the plink-ng submodule:
#   1. plink_libdeflate  — libdeflate (C), used by pgenlib for BGZF support
#   2. plink_zstd        — vendored zstd (C), used by pgenlib for .zst files
#   3. pgenlib           — plink2 genomics library (C++17)
#
# Hidden visibility is set on all three libraries to prevent symbol collisions
# with DuckDB's own bundled zstd. Without this, duplicate zstd symbols would
# cause undefined behavior at runtime.
# ===========================================================================

# --- Paths ---
set(PLINK_NG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/plink-ng/2.0)
set(PGENLIB_DIR ${PLINK_NG_DIR}/include)
set(LIBDEFLATE_DIR ${PLINK_NG_DIR}/libdeflate)
set(ZSTD_DIR ${PLINK_NG_DIR}/zstd)
set(SIMDE_DIR ${PLINK_NG_DIR}/simde)

# Fail early if submodule is not initialized
if(NOT EXISTS ${PGENLIB_DIR}/pgenlib_read.cc)
    message(FATAL_ERROR
        "plink-ng submodule not initialized. Run:\n"
        "  git submodule update --init third_party/plink-ng")
endif()

# --- libdeflate (C library) ---
add_library(plink_libdeflate STATIC
    ${LIBDEFLATE_DIR}/lib/adler32.c
    ${LIBDEFLATE_DIR}/lib/crc32.c
    ${LIBDEFLATE_DIR}/lib/deflate_compress.c
    ${LIBDEFLATE_DIR}/lib/deflate_decompress.c
    ${LIBDEFLATE_DIR}/lib/gzip_compress.c
    ${LIBDEFLATE_DIR}/lib/gzip_decompress.c
    ${LIBDEFLATE_DIR}/lib/utils.c
    ${LIBDEFLATE_DIR}/lib/zlib_compress.c
    ${LIBDEFLATE_DIR}/lib/zlib_decompress.c
    ${LIBDEFLATE_DIR}/lib/arm/arm_cpu_features.c
    ${LIBDEFLATE_DIR}/lib/x86/x86_cpu_features.c
)
target_include_directories(plink_libdeflate PUBLIC ${LIBDEFLATE_DIR})
target_compile_definitions(plink_libdeflate PRIVATE LIBDEFLATE_STATIC)
set_target_properties(plink_libdeflate PROPERTIES
    C_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON
)

# --- zstd (C library, pgenlib's vendored copy) ---
# Hidden visibility is critical to avoid symbol collisions with DuckDB's bundled zstd.
# GLOB is acceptable here because the submodule is pinned; re-run cmake if pin changes.
file(GLOB PLINK_ZSTD_COMMON_SRCS ${ZSTD_DIR}/lib/common/*.c)
file(GLOB PLINK_ZSTD_DECOMPRESS_SRCS ${ZSTD_DIR}/lib/decompress/*.c)
file(GLOB PLINK_ZSTD_COMPRESS_SRCS ${ZSTD_DIR}/lib/compress/*.c)
add_library(plink_zstd STATIC
    ${PLINK_ZSTD_COMMON_SRCS}
    ${PLINK_ZSTD_DECOMPRESS_SRCS}
    ${PLINK_ZSTD_COMPRESS_SRCS}
)
# PRIVATE: don't leak vendored zstd headers to extension targets (collision risk
# with DuckDB's own bundled zstd headers).
target_include_directories(plink_zstd PRIVATE
    ${ZSTD_DIR}/lib
    ${ZSTD_DIR}/lib/common
)
target_compile_definitions(plink_zstd PRIVATE ZSTD_DISABLE_ASM)
set_target_properties(plink_zstd PROPERTIES
    C_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON
)

# --- musl libc compatibility ---
# rawmemchr is a GNU extension that plink-ng uses when _GNU_SOURCE is defined.
# glibc provides it, but musl (which also defines _GNU_SOURCE) does not.
# Detect at configure time and force-include a shim header if missing.
include(CheckSymbolExists)
set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
check_symbol_exists(rawmemchr "string.h" HAVE_RAWMEMCHR)

# --- pgenlib (C++ library) ---
# Read-path only — pgenlib_write.cc, plink2_bitmap.cc, plink2_fmath.cc,
# plink2_stats.cc are intentionally excluded.
add_library(pgenlib STATIC
    ${PGENLIB_DIR}/plink2_base.cc
    ${PGENLIB_DIR}/plink2_bits.cc
    ${PGENLIB_DIR}/plink2_bgzf.cc
    ${PGENLIB_DIR}/plink2_htable.cc
    ${PGENLIB_DIR}/plink2_memory.cc
    ${PGENLIB_DIR}/plink2_string.cc
    ${PGENLIB_DIR}/plink2_text.cc
    ${PGENLIB_DIR}/plink2_thread.cc
    ${PGENLIB_DIR}/plink2_zstfile.cc
    ${PGENLIB_DIR}/pgenlib_misc.cc
    ${PGENLIB_DIR}/pgenlib_read.cc
    ${PGENLIB_DIR}/pgenlib_ffi_support.cc
    ${PGENLIB_DIR}/pvar_ffi_support.cc
)
target_include_directories(pgenlib PUBLIC
    ${PGENLIB_DIR}
    ${LIBDEFLATE_DIR}
    ${SIMDE_DIR}
)
target_include_directories(pgenlib PRIVATE
    ${ZSTD_DIR}/lib
    ${ZSTD_DIR}/lib/common
)
target_link_libraries(pgenlib PRIVATE plink_libdeflate plink_zstd)
target_compile_definitions(pgenlib PRIVATE
    ZSTD_DISABLE_ASM
    LIBDEFLATE_STATIC
)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries(pgenlib PRIVATE ZLIB::ZLIB Threads::Threads)
set_target_properties(pgenlib PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
)
if(NOT HAVE_RAWMEMCHR)
    target_compile_options(pgenlib PRIVATE
        -include ${CMAKE_CURRENT_SOURCE_DIR}/src/musl_compat.h
    )
endif()

# --- Link pgenlib and its dependencies to extension targets ---
# Deps are listed explicitly because pgenlib links them PRIVATE (won't propagate
# in static builds) and DuckDB's macros require the plain target_link_libraries
# signature (no PRIVATE/PUBLIC keyword).
target_link_libraries(${EXTENSION_NAME} pgenlib plink_libdeflate plink_zstd ZLIB::ZLIB Threads::Threads)
target_link_libraries(${LOADABLE_EXTENSION_NAME} pgenlib plink_libdeflate plink_zstd ZLIB::ZLIB Threads::Threads)
target_include_directories(${EXTENSION_NAME} PRIVATE ${PGENLIB_DIR})
target_include_directories(${LOADABLE_EXTENSION_NAME} PRIVATE ${PGENLIB_DIR})

install(
  TARGETS ${EXTENSION_NAME} pgenlib plink_libdeflate plink_zstd
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
