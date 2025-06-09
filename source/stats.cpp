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

#include <numbers>

#include "stats.hpp"

#include "logger.hpp"

namespace gamba
{

double stats::update_cputime{0.0};
double stats::update_walltime{0.0};

double stats::select_cputime{0.0};
double stats::select_walltime{0.0};

double stats::matrix_cputime{0.0};
double stats::matrix_walltime{0.0};

double stats::symbolic_cputime{0.0};
double stats::symbolic_walltime{0.0};

double stats::convert_cputime{0.0};
double stats::convert_walltime{0.0};

double stats::linalg_cputime{0.0};
double stats::linalg_walltime{0.0};

double stats::insert_cputime{0.0};
double stats::insert_walltime{0.0};

double stats::reduce_cputime{0.0};
double stats::reduce_walltime{0.0};

double stats::linalg_interred_cputime{0.0};
double stats::linalg_interred_walltime{0.0};

double stats::reconstruct_cputime{0.0};
double stats::reconstruct_walltime{0.0};

double stats::overall_cputime{0.0};
double stats::overall_walltime{0.0};

ssize_t stats::spairs_reduced{0};

ssize_t stats::gm_criteria{0};

ssize_t stats::redundant_elem{0};

size_t stats::rows_reduced{0};

size_t stats::zero_reductions{0};

size_t stats::max_size_bht{0};
size_t stats::max_size_sht{0};
size_t stats::max_size_mht{0};

size_t stats::num_primes{0};

void stats::print_timings()
{
    log::print(log::INFO1, "\n┌{0:─^{1}}┐\n", " TIMINGS ", 38);

    log::print(log::INFO1, "│ overall (wall) {:>17.2f} sec │\n",
               overall_walltime);

    log::print(log::INFO1, "│ overall (cpu) {:>11.2f} sec {:5.1f}x │\n",
               overall_cputime, overall_cputime / overall_walltime);

    log::print(log::INFO1, "│ matrix const. {:>11.2f} sec {:5.1f}% │\n",
               matrix_walltime + select_walltime,
               100.0 * (matrix_walltime + select_walltime) / overall_walltime);

    log::print(log::INFO1, "│ symbolic prep. {:>10.2f} sec {:5.1f}% │\n",
               symbolic_walltime, 100.0 * symbolic_walltime / overall_walltime);

    log::print(log::INFO1, "│ convert cols. {:>11.2f} sec {:5.1f}% │\n",
               convert_walltime, 100.0 * convert_walltime / overall_walltime);

    log::print(log::INFO1, "│ linear algebra {:10.2f} sec {:5.1f}% │\n",
               linalg_walltime, 100.0 * linalg_walltime / overall_walltime);

    log::print(log::INFO1, "│ insert rows {:13.2f} sec {:5.1f}% │\n",
               insert_walltime, 100.0 * insert_walltime / overall_walltime);

    log::print(log::INFO1, "│ update spairs {:11.2f} sec {:5.1f}% │\n",
               update_walltime, 100.0 * update_walltime / overall_walltime);

    log::print(log::INFO1, "│ reduce basis {:12.2f} sec {:5.1f}% │\n",
               reduce_walltime, 100.0 * reduce_walltime / overall_walltime);

    if (reconstruct_walltime > 0.0)
    {
        log::print(log::INFO1, "│ rational recon. {:9.2f} sec {:5.1f}% │\n",
                   reconstruct_walltime,
                   100.0 * reconstruct_walltime / overall_walltime);
    }

    log::print(log::INFO1, "└{0:─^{1}}┘\n", "", 38);
}

void stats::print_statistics()
{
    log::print(log::INFO1, "\n┌{0:─^{1}}┐\n", " STATISTICS ", 38);

    log::print(log::INFO1, "│ num. spairs reduced {:>16} │\n", spairs_reduced);

    log::print(log::INFO1, "│ num. GM criterion {:>18} │\n", gm_criteria);

    log::print(log::INFO1, "│ num. redundant elem. {:>15} │\n", redundant_elem);

    log::print(log::INFO1, "│ num. rows reduced {:>18} │\n", rows_reduced);

    log::print(log::INFO1, "│ num. zero reductions {:>15} │\n",
               zero_reductions);

    log::print(log::INFO1, "│ max. size basis ht {0:>15}{1:} │\n", "2^",
               std::ceil(std::log(static_cast<double>(max_size_bht))
                         / std::numbers::ln2));

    log::print(log::INFO1, "│ max. size spair ht {0:>15}{1:} │\n", "2^",
               std::ceil(std::log(static_cast<double>(max_size_sht))
                         / std::numbers::ln2));

    log::print(log::INFO1, "│ max. size matrix ht {0:>14}{1:} │\n", "2^",
               std::ceil(std::log(static_cast<double>(max_size_mht))
                         / std::numbers::ln2));

    if (num_primes > 0)
    {
        log::print(log::INFO1, "│ num. primes used {0:>19} │\n", num_primes);
    }

    log::print(log::INFO1, "└{0:─^{1}}┘\n", "", 38);
}

void stats::reset_timings()
{
    update_cputime  = 0.0;
    update_walltime = 0.0;

    select_cputime  = 0.0;
    select_walltime = 0.0;

    matrix_cputime  = 0.0;
    matrix_walltime = 0.0;

    symbolic_cputime = 0.0;

    symbolic_walltime = 0.0;

    convert_cputime  = 0.0;
    convert_walltime = 0.0;

    linalg_cputime  = 0.0;
    linalg_walltime = 0.0;

    insert_cputime  = 0.0;
    insert_walltime = 0.0;

    reduce_cputime  = 0.0;
    reduce_walltime = 0.0;

    linalg_interred_cputime  = 0.0;
    linalg_interred_walltime = 0.0;

    reconstruct_cputime  = 0.0;
    reconstruct_walltime = 0.0;

    overall_cputime  = 0.0;
    overall_walltime = 0.0;
}

void stats::reset_statistics()
{
    spairs_reduced = 0;

    gm_criteria = 0;

    redundant_elem = 0;

    rows_reduced = 0;

    zero_reductions = 0;

    max_size_bht = 0;
    max_size_sht = 0;
    max_size_mht = 0;

    num_primes = 0;
}

}  // namespace gamba
