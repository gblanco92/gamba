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

#include <chrono>

#include "basis.hpp"
#include "linalg_v3.hpp"
#include "matrix.hpp"
#include "reduce.hpp"
#include "update.hpp"

namespace gamba
{

namespace f4
{

void print_column_names();

void print_bottom_line();

}  // namespace f4

template <class MonomialOrder, class CoefficientType>
void f4_main(polynomial_basis<CoefficientType>& basis)
{
    using coefficient_type = CoefficientType;
    using monomial_order   = MonomialOrder;
    using linalg_type      = linalg_v3<coefficient_type>;

    spair_set spairs;

    /* generate first spairs and update redundant basis elements */
    update_f4<monomial_order>(spairs, basis, 0);

    /* reuse allocated memory in matrix across rounds */
    matrix_f4 matrix{basis};

    /* reuse allocated memory in linear algebra across rounds */
    linalg_type echelon_engine3{basis.field};

    f4::print_column_names();

    /* main f4 loop */
    while (not spairs.queue.empty())
    {
        /* timings */
        auto const start_walltime = std::chrono::system_clock::now();

        size_t const prev_num_gens = basis.num_gens();

        auto const [num_spairs, round_degree] =
            select_spairs<monomial_order>(spairs);

        // print_memory_usage(spairs.memory_usage());

        matrix.insert_spairs(spairs, num_spairs, basis);

        matrix.symbolic_preprocessing(basis);

        matrix.convert_monomials_to_columns(monomial_order{});

        // print_memory_usage(matrix.memory_usage());

        echelon_engine3.initialize(matrix);

        auto const echl_time = echelon_engine3.reduce();

        // print_memory_usage(echelon_engine3.memory_usage());

        print_time(echl_time);

        echelon_engine3.extract_new_rows(matrix);

        basis.insert_new_rows_echelon(matrix, monomial_order{});

        // print_memory_usage(basis.memory_usage());

        update_f4<monomial_order>(spairs, basis, prev_num_gens);

        /* timings */
        auto const end_walltime = std::chrono::system_clock::now();

        print_memory_usage();

        print_time(end_walltime - start_walltime, /* new_line = */ true);

        if (basis.is_trivial())
            break;
    }

    f4::print_bottom_line();

    basis.remove_redundant_gens();

    /* matrix has been cleared by basis */
    if (not params::no_reduce)
        reduce<monomial_order>(basis, matrix);
}

}  // namespace gamba
