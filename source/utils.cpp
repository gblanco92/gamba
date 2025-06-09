/*   GamBa: a Groebner basis engine
 *   Copyright (C) 2024 Guillem Blanco
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>. */

#include "utils.hpp"

#include <gmp.h>

#include "thirdparty/getRSS.hpp"

namespace gamba
{

std::string compiler_version_string()
{
    constexpr std::string_view const cxx =
#ifdef __clang__
        "clang++";
#else
        "g++";
#endif

    constexpr int32_t const cxx_ver_major =
#ifdef __clang__
        __clang_major__;
#else
        __GNUC__;
#endif

    constexpr int32_t const cxx_ver_minor =
#ifdef __clang__
        __clang_minor__;
#else
        __GNUC_MINOR__;
#endif

    constexpr int32_t const cxx_ver_patch =
#ifdef __clang__
        __clang_patchlevel__;
#else
        __GNUC_PATCHLEVEL__;
#endif

    return fmt::format("{}-{}.{}.{}", cxx, cxx_ver_major, cxx_ver_minor,
                       cxx_ver_patch);
}

std::string libgmp_version_string()
{
    return fmt::format("gmp-{}.{}.{}", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR,
                       __GNU_MP_VERSION_PATCHLEVEL);
}

std::string libflint_version_string()
{
    return fmt::format("flint-{}.{}.{}", __FLINT_VERSION, __FLINT_VERSION_MINOR,
                       __FLINT_VERSION_PATCHLEVEL);
}

void print_memory_usage(double const mem_usage)
{
    if (mem_usage < 1024.0)
        log::print(log::INFO2, "{:10.2f} MiB", mem_usage);
    else
        log::print(log::INFO2, "{:10.2f} GiB", mem_usage / 1024.0);
}

void print_memory_usage()
{
    /* get peak RSS memory usage in mebibytes */
    double const mem_usage =
        static_cast<double>(getPeakRSS()) / 1024.0 / 1024.0;

    print_memory_usage(mem_usage);
}

void print_time(std::chrono::duration<double> const time, bool const new_line)
{
    log::print(log::INFO2, "{:10.2f} sec", time.count());

    if (new_line)
        log::print(log::INFO2, " │\n");

    ::fflush(stdout);
}

}  // namespace gamba
