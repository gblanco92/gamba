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

#pragma once

#include "basis.hpp"
#include "f4.hpp"
#include "modular.hpp"
#include "monomial.hpp"
#include "stats.hpp"

namespace gamba
{

template <class CoefficientType>
generators_data groebner_basis_main(generators_data const& input_data,
                                    monomial_order const& mon_order)
{
    using coeff_type = CoefficientType;
    using basis_type = polynomial_basis<coeff_type>;

    /* set monomial_base static data */
    monomial_base::exp_size = mon_order.exponent_size(input_data.num_vars);

    if (mon_order.is_block_order())
        monomial_base::block_size = params::num_elim_vars + 1;

    monomial_base::initialize_weights();

    /* initialize basis from input generators */
    basis_type basis{input_data.num_vars, input_data.field_char};
    basis.import_generators(input_data, mon_order);

    GAMBA_DEBUG(basis.print_generators());

    basis.print_info();

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    constexpr bool const rational_coeffs =
        std::is_same_v<coeff_type, fmpq_class>;

    if constexpr (rational_coeffs)
    {
        f4_modular(basis, mon_order);
    }
    else
    {
        f4_main(basis, mon_order, nullptr);
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::overall_walltime =
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::overall_cputime =
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    stats::print_timings();
    stats::print_statistics();

    generators_data output_data;
    basis.export_generators(output_data, params::lead_mons);

    /* copy variable names from input data */
    output_data.var_names = input_data.var_names;

    return output_data;
}

}  // namespace gamba
