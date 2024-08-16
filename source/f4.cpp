/*   GamBa: a Groebner basis engine
 *   Copyright (C) 2023 Guillem Blanco
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

#include <format>
#include <iostream>

#include "getRSS.hpp"

namespace gamba
{

void print_column_names()
{
    std::cout << std::endl;
    std::cout << std::format("{}{:>12}{:>18}{:>15}{:>20}{:>19}{:>14}{:>14}",
                             "deg", "spairs", "matrix", "density", "new gens",
                             "ech. time", "mem. usage", "total time")
              << std::endl;
    std::cout << "---------------------------------------------------------";
    std::cout << "----------------------------------------------------------"
              << std::endl;
}

void print_memory_usage(double const mem_usage)
{
    if (mem_usage < 1024.0)
        std::cout << std::format("{:10.2f} MiB", mem_usage);
    else
        std::cout << std::format("{:10.2f} GiB", mem_usage / 1024.0);
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
    std::cout << std::format("{:10.2f} sec", time.count()) << std::flush;

    if (new_line)
        std::cout << std::endl;
}

void print_bottom_line()
{
    std::cout << "=========================================================";
    std::cout << "=========================================================="
              << std::endl;
}

}  // namespace gamba
