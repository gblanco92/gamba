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

#pragma once

namespace gamba
{

struct stats
{
    static void init();

    static void reset_statistics();

    static void reset_timings();

    static void print_statistics();

    static void print_timings();

    /* GamBa timings */
    static double update_cputime;
    static double update_walltime;

    static double select_cputime;
    static double select_walltime;

    static double matrix_cputime;
    static double matrix_walltime;

    static double symbolic_cputime;
    static double symbolic_walltime;

    static double convert_cputime;
    static double convert_walltime;

    static double linalg_cputime;
    static double linalg_walltime;

    static double insert_cputime;
    static double insert_walltime;

    static double reduce_cputime;
    static double reduce_walltime;

    static double linalg_interred_cputime;
    static double linalg_interred_walltime;

    static double reconstruct_cputime;
    static double reconstruct_walltime;

    static double overall_cputime;
    static double overall_walltime;

    /* GamBa stats */
    static ssize_t spairs_reduced;

    static ssize_t gm_criteria;

    static ssize_t redundant_elem;

    static size_t rows_reduced;

    static size_t zero_reductions;

    static size_t max_size_bht;
    static size_t max_size_sht;
    static size_t max_size_mht;

    static size_t num_primes;
};

}  // namespace gamba
