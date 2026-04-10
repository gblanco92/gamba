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
#include "linalg.hpp"
#include "utils.hpp"

namespace gamba
{

template <class CoefficientType>
void reduce(polynomial_basis<CoefficientType>& basis,
            monomial_order const& mon_order,
            matrix_f4& matrix)
{
    using coeff_type = CoefficientType;

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    log::print(log::INFO2, "│ reduce basis       ");

    matrix.insert_generators_reduce(basis);

    matrix.symbolic_preprocessing(basis);
    /* do not count the reduced rows from the reduction phase */
    stats::rows_reduced -= matrix.num_bottom_rows();

    matrix.convert_monomials_to_columns(mon_order);

    linalg<coeff_type> interreduce_engine{basis.field, matrix};
    interreduce_engine.interreduce(matrix);

    basis.clear();
    basis.insert_new_rows_reduce(matrix, mon_order);

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::reduce_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::reduce_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    print_memory_usage(log::INFO2);

    std::chrono::duration<double> const time = end_walltime - start_walltime;

    log::print(log::INFO2, "{:10.2f} sec │\n", time.count());
    log::print(log::INFO2, "└{0:─^{1}}┘\n", "", 118);
}

}  // namespace gamba
