cmake_minimum_required(VERSION 3.20)

project(
    GamBa
    VERSION 0.2.0
    DESCRIPTION "A Groebner basis engine"
    HOMEPAGE_URL "https://github.com/gblanco92/GamBa"
    LANGUAGES CXX
)

# Includes
include(cmake/dev-mode.cmake)
include(cmake/prelude.cmake)
include(FetchContent)

set(pre_configure_dir "${PROJECT_SOURCE_DIR}/source")
set(post_configure_dir "${PROJECT_SOURCE_DIR}/source")
include(cmake/CheckGit.cmake)
CheckGitSetup()

# Include cmake/ into CMake path
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)

if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# ---- Dependencies ----

FetchContent_Declare(
        CpuInfo
        QUIET
        GIT_REPOSITORY  "https://github.com/pytorch/cpuinfo/"
        GIT_TAG         "main"
        BINARY_DIR      ${CMAKE_BINARY_DIR}/cpuinfo/
)

set(CPUINFO_BUILD_TOOLS OFF CACHE BOOL "Do not build command-line tools")
set(CPUINFO_BUILD_UNIT_TESTS OFF CACHE BOOL "Do not build cpuinfo unit tests")
set(CPUINFO_BUILD_MOCK_TESTS OFF CACHE BOOL "Do not build cpuinfo mock tests")
set(CPUINFO_BUILD_BENCHMARKS OFF CACHE BOOL "Do not build cpuinfo micro-benchmarks")
set(CPUINFO_LIBRARY_TYPE "static" CACHE BOOL "Built only an static library")

FetchContent_MakeAvailable(CpuInfo)

#target_compile_features(CpuInfo PRIVATE cxx_std_23)
#set_target_properties(CpuInfo PROPERTIES CXX_STANDARD_REQUIRED ON)
#set_target_properties(CpuInfo PROPERTIES CXX_EXTENSIONS OFF)

FetchContent_Declare(
        cli11
        QUIET
        GIT_REPOSITORY  "https://github.com/CLIUtils/CLI11"
        GIT_TAG         v2.5.0
)

FetchContent_MakeAvailable(cli11)

FetchContent_Declare(
        fmt
        GIT_REPOSITORY  "https://github.com/fmtlib/fmt"
        GIT_TAG         11.1.4
        BINARY_DIR      ${CMAKE_BINARY_DIR}/fmt/
) 

FetchContent_MakeAvailable(fmt)

# ----------------------

set(GamBa_SRCS source/basis.cpp
               source/base_basis.cpp
               source/f4.cpp
               source/gamba.cpp
               source/gb_main.cpp
               source/thirdparty/getRSS.cpp
               source/io.cpp
               source/kernel/avx2/u16/echelon_u16.cpp
               source/kernel/avx2/u32/echelon_u32.cpp
               source/kernel/avx2/find.cpp
               source/kernel/avx2/u16/misc_u16.cpp
               source/kernel/avx2/u32/misc_u32.cpp
               source/kernel/magic.cpp
               source/kernel/avx2/u16/modular_u16.cpp
               source/kernel/avx2/u32/modular_u32.cpp
               source/kernel/avx2/u16/reduce_u16.cpp
               source/kernel/avx2/u32/reduce_u32.cpp
               source/kernel/avx2/saxpy.cpp
               source/linalg.cpp
               source/linalg_v3.cpp
               source/linalg_pack.cpp
               source/main.cpp
               source/matrix.cpp
               source/memalign.cpp
               source/modular.cpp
               source/params.cpp
               source/profile.cpp
               source/reduce.cpp
               source/spair.cpp
               source/stats.cpp
               source/update.cpp
               source/utils.cpp)

# ---- Declare library ----

##add_library(GamBa_lib OBJECT source/gamba.cpp)
##target_include_directories(GamBa_lib PUBLIC ${PROJECT_SOURCE_DIR}/source)
##set_target_properties(GamBa_lib PROPERTIES OUTPUT_NAME gamba)

# ---- Declare executable ----

add_executable(GamBa_exe ${GamBa_SRCS})
target_include_directories(GamBa_exe PUBLIC ${PROJECT_SOURCE_DIR}/source)
set_target_properties(GamBa_exe PROPERTIES OUTPUT_NAME gamba)

# ---- Link libraries ----

target_link_libraries(GamBa_exe PUBLIC git_version)
target_link_libraries(GamBa_exe PRIVATE cpuinfo CLI11::CLI11 fmt::fmt)
target_link_libraries(GamBa_exe PRIVATE libgmp.a libgmpxx.a libmpfr.a libflint.a)

if (CMAKE_BUILD_TYPE STREQUAL "Release" AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_link_libraries(GamBa_exe PRIVATE libc++.a libc++abi.a)
endif()

#target_link_libraries(GamBa_exe PRIVATE tcmalloc_minimal.a)
#target_link_libraries(GamBa_exe PRIVATE jemalloc.a)

# ---- Compile features ----

target_compile_features(GamBa_exe PRIVATE cxx_std_23)
set_target_properties(GamBa_exe PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(GamBa_exe PROPERTIES CXX_EXTENSIONS OFF)

target_compile_features(fmt PRIVATE cxx_std_23)
set_target_properties(fmt PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(fmt PROPERTIES CXX_EXTENSIONS OFF)

if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "ReleaseDebug")
  set_target_properties(GamBa_exe PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
  set_target_properties(fmt PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)

  set_target_properties(GamBa_exe PROPERTIES UNITY_BUILD ON)
  set_target_properties(GamBa_exe PROPERTIES UNITY_BUILD_MODE GROUP)
  set_target_properties(GamBa_exe PROPERTIES UNITY_BUILD_UNIQUE_ID "GAMBA_UNITY_BUILD")
  set_source_files_properties(${GamBa_SRCS} PROPERTIES UNITY_GROUP "gamba")
endif()

# ---- Mandatory flags ----

target_compile_options(GamBa_exe PRIVATE -march=core-avx2 -mno-vzeroupper)

# ---- Precompiled headers ----

target_precompile_headers(GamBa_exe PRIVATE
                          <algorithm>
                          <bit>
                          <cassert>
                          <chrono>
                          <CLI/CLI.hpp>
                          <climits>
                          <cmath>
                          <cpuinfo.h>
                          <concepts>
                          <cstddef>
                          <cstdint>
                          <cstdio>
                          <cstdlib>
                          <cstring>
                          <ctime>
                          <fcntl.h>
                          <flint/flint.h>
                          <flint/fmpz.h>
                          <flint/fmpq.h>
                          <fmt/base.h>
                          <fmt/format.h>
                          <fmt/ranges.h>
                          <fmt/chrono.h>
                          <fmt/std.h>
                          <fmt/compile.h>
                          <fmt/color.h>
                          <fmt/os.h>
                          <fmt/ostream.h>
                          <fmt/args.h>
                          <fmt/printf.h>
                          <fmt/xchar.h>
                          <format>
                          <fstream>
                          <gmpxx.h>
                          <immintrin.h>
                          <iomanip>
                          <iostream>
                          <limits>
                          <memory>
                          <new>
                          <numeric>
                          <random>
                          <ranges>
                          <string>
                          <sys/mman.h>
                          <sys/resource.h>
                          <sys/stat.h>
                          <sys/types.h>
                          <type_traits>
                          <unistd.h>
                          <utility>
                          <vector>
                          <xmmintrin.h>)