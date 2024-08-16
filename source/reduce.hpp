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

#include "basis.hpp"
#include "f4.hpp"
#include "linalg.hpp"
#include "params.hpp"

namespace gamba
{

template <class CoefficientType, class MonomialOrder>
void reduce(polynomial_basis<CoefficientType, MonomialOrder>& basis,
            matrix_f4<CoefficientType, MonomialOrder>& matrix,
            gamba_params const& params)
{
    using coefficient_type = CoefficientType;
    using monomial_order   = MonomialOrder;
    using matrix_type      = matrix_f4<coefficient_type, monomial_order>;

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    if (params.rounds_info)
        std::cout << "reduce basis      ";

    matrix.insert_generators_reduce(basis);

    matrix.symbolic_preprocessing(basis);
    /* do not count the reduced rows from the reduction phase */
    stats.rows_reduced -= matrix.num_bottom_rows();

    matrix.convert_monomials_to_columns_reduce(params);

    linalg<matrix_type> reduce_engine{matrix};
    reduce_engine.reduce(matrix, params);

    basis.clear();
    basis.insert_new_rows_reduce(matrix);

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.reduce_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.reduce_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    if (params.rounds_info)
    {
        std::cout << std::format("{:10.2f} sec", stats.linalg_reduce_walltime)
                  << std::flush;

        print_memory_usage();

        print_time(end_walltime - start_walltime);

        std::cout << std::endl;
    }
}

}  // namespace gamba
