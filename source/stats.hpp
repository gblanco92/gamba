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

#include <sys/types.h>

namespace gamba
{

struct f4_statistics
{
    /* timings */
    double update_cputime{0.0};
    double update_walltime{0.0};

    double select_cputime{0.0};
    double select_walltime{0.0};

    double matrix_cputime{0.0};
    double matrix_walltime{0.0};

    double symbolic_cputime{0.0};
    double symbolic_walltime{0.0};

    double convert_cputime{0.0};
    double convert_walltime{0.0};

    double linalg_cputime{0.0};
    double linalg_walltime{0.0};

    double insert_cputime{0.0};
    double insert_walltime{0.0};

    double reduce_cputime{0.0};
    double reduce_walltime{0.0};

    double linalg_reduce_cputime{0.0};
    double linalg_reduce_walltime{0.0};

    double overall_cputime{0.0};
    double overall_walltime{0.0};

    /* F4 data */
    ssize_t spairs_reduced{0};

    ssize_t gm_criteria{0};

    ssize_t redundant_elem{0};

    size_t rows_reduced{0};

    size_t zero_reductions{0};

    size_t max_size_bht{0};
    size_t max_size_sht{0};
    size_t max_size_mht{0};
};

extern f4_statistics stats;

void print_statistics();

void print_timings();

}  // namespace gamba
