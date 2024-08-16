cmake_minimum_required(VERSION 3.20)

project(
    GamBa
    VERSION 0.1.0
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

FetchContent_Declare(
        cli11_proj
        QUIET
        GIT_REPOSITORY  "https://github.com/CLIUtils/CLI11"
        GIT_TAG         v2.4.2
)

FetchContent_MakeAvailable(cli11_proj)

# ---- Declare library ----

##add_library(
##    GamBa_lib OBJECT
##    source/gamba.cpp
##)
##
##target_include_directories(
##    GamBa_lib ${warning_guard}
##    PUBLIC
##    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
##)
##
##target_compile_features(GamBa_lib PUBLIC cxx_std_23)

# ---- Declare executable ----

add_executable(GamBa_exe source/debug.cpp
                         source/f4.cpp
                         source/gamba.cpp
                         source/getRSS.cpp
                         source/io.cpp
                         source/kernel/echelon_u16.cpp
                         source/kernel/echelon_u32.cpp
                         source/kernel/find.cpp
                         source/kernel/misc_u16.cpp
                         source/kernel/misc_u32.cpp
                         source/kernel/modular_u16.cpp
                         source/kernel/modular_u32.cpp
                         source/kernel/reduce_u16.cpp
                         source/kernel/reduce_u32.cpp
                         source/kernel/saxpy.cpp
                         source/main.cpp
                         source/params.cpp
                         source/profile.cpp
                         source/stats.cpp)

target_include_directories(GamBa_exe PUBLIC ${PROJECT_SOURCE_DIR}/source)

add_executable(GamBa::exe ALIAS GamBa_exe)

# ---- Target properties ----

set_target_properties(GamBa_exe PROPERTIES OUTPUT_NAME gamba)

target_link_libraries(GamBa_exe PUBLIC git_version)

target_link_libraries(GamBa_exe PRIVATE cpuinfo CLI11::CLI11)

target_link_libraries(GamBa_exe PRIVATE libgmp.a libgmpxx.a)

if (CMAKE_BUILD_TYPE STREQUAL "Release" AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_link_libraries(GamBa_exe PRIVATE libc++.a libc++abi.a)
endif()

#target_link_libraries(GamBa_exe PRIVATE tcmalloc_minimal.a)

#target_link_libraries(GamBa_exe PRIVATE jemalloc.a)

target_compile_features(GamBa_exe PRIVATE cxx_std_23)

set_target_properties(GamBa_exe PROPERTIES CXX_STANDARD_REQUIRED ON)

set_target_properties(GamBa_exe PROPERTIES CXX_EXTENSIONS OFF)

if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "ReleaseDebug")
  set_target_properties(GamBa_exe PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
endif()

# ---- Target flags ----

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=core-avx2 -mno-vzeroupper")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

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
                          <cpuinfo.h>
                          <cstddef>
                          <cstdint>
                          <cstdio>
                          <cstdlib>
                          <cstring>
                          <ctime>
                          <fcntl.h>
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