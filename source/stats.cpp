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

#include "stats.hpp"

namespace gamba
{

/* use a global singleton for logging statistics during the computation,
 * not the best approach but it will work for now */
f4_statistics stats{};

void print_timings()
{
    std::cout << std::endl;
    std::cout << "*************** TIMINGS ***************" << std::endl;
    std::cout << std::format("overall (wall) {:>20.2f} sec",
                             stats.overall_walltime)
              << std::endl;

    std::cout << std::format("overall (cpu) {:>14.2f} sec {:5.1f}x",
                             stats.overall_cputime,
                             stats.overall_cputime / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "matrix const. {:>14.2f} sec {:5.1f}%",
        stats.matrix_walltime + stats.select_walltime,
        100.0 * (stats.matrix_walltime + stats.select_walltime)
            / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "symbolic prep. {:>13.2f} sec {:5.1f}%", stats.symbolic_walltime,
        100.0 * stats.symbolic_walltime / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "convert cols {:15.2f} sec {:5.1f}%", stats.convert_walltime,
        100.0 * stats.convert_walltime / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "linear algebra {:13.2f} sec {:5.1f}%", stats.linalg_walltime,
        100.0 * stats.linalg_walltime / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "insert rows {:16.2f} sec {:5.1f}%", stats.insert_walltime,
        100.0 * stats.insert_walltime / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "update spairs {:14.2f} sec {:5.1f}%", stats.update_walltime,
        100.0 * stats.update_walltime / stats.overall_walltime)
              << std::endl;

    std::cout << std::format(
        "reduce basis {:15.2f} sec {:5.1f}%", stats.reduce_walltime,
        100.0 * stats.reduce_walltime / stats.overall_walltime)
              << std::endl;

    std::cout << "***************************************" << std::endl;
}

void print_statistics()
{
    std::cout << std::endl;
    std::cout << "*************** F4 DATA ***************" << std::endl;

    std::cout << std::format("num. spairs reduced {:>19}", stats.spairs_reduced)
              << std::endl;

    std::cout << std::format("num. GM criterion {:>21}", stats.gm_criteria)
              << std::endl;

    std::cout << std::format("num. redundant elem. {:>18}",
                             stats.redundant_elem)
              << std::endl;

    std::cout << std::format("num. rows reduced {:>21}", stats.rows_reduced)
              << std::endl;

    std::cout << std::format("num. zero reductions {:>18}",
                             stats.zero_reductions)
              << std::endl;

    std::cout << std::format(
        "max. size basis ht {0:>18}{1:}", "2^",
        std::ceil(std::log(static_cast<double>(stats.max_size_bht))
                  / std::log(2)))
              << std::endl;

    std::cout << std::format(
        "max. size spair ht {0:>18}{1:}", "2^",
        std::ceil(std::log(static_cast<double>(stats.max_size_sht))
                  / std::log(2)))
              << std::endl;

    std::cout << std::format(
        "max. size matrix ht {0:>17}{1:}", "2^",
        std::ceil(std::log(static_cast<double>(stats.max_size_mht))
                  / std::log(2)))
              << std::endl;

    std::cout << "***************************************" << std::endl;
}

}  // namespace gamba
