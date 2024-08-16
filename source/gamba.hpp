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

#include <limits>
#include <stdexcept>

#include <gmpxx.h>

#include "basis.hpp"
#include "debug.hpp"
#include "f4.hpp"
#include "monomial.hpp"
#include "params.hpp"
#include "stats.hpp"

namespace gamba
{

generators_data groebner_basis(generators_data const& data,
                               gamba_params& params);

template <class MonomialOrder>
generators_data groebner_basis_coeff_type(generators_data const& data,
                                          gamba_params& params)
{
    if (data.field_char == 2)
    {
        throw std::runtime_error("Field characteristic = 2 not supported.");
    }

    if (data.field_char == 0)
    {
        // TODO: not implemented yet
        return generators_data{};
    }

    if (data.field_char < std::numeric_limits<uint8_t>::max())
    {
        return gb_main<uint8_t, MonomialOrder>(data, params);
    }

    if (data.field_char < std::numeric_limits<uint16_t>::max())
    {
        try
        {
            return gb_main<uint16_t, MonomialOrder>(data, params);
        }
        catch (std::overflow_error const&)
        {
            if (params.verbose > 0)
            {
                std::cout << " Index overflow detected. Restarting..."
                          << std::endl;
            }

            params.basis_info = false;

            return gb_main<uint32_t, MonomialOrder>(data, params);
        }
    }

    /* max. allowed field characteristic is 2^31 - 1 */
    if (data.field_char <= std::numeric_limits<int32_t>::max())
    {
        return gb_main<uint32_t, MonomialOrder>(data, params);
    }

    throw std::runtime_error("Field characteristic too large.");
}

template <class CoefficientType, class MonomialOrder>
generators_data gb_main(generators_data const& input_data,
                        gamba_params const& params)
{
    using monomial_order = MonomialOrder;
    using coeff_type     = CoefficientType;
    using basis_type     = polynomial_basis<coeff_type, monomial_order>;

    /* set monomial_order & monomial_base static data */
    if constexpr (is_block_order_v<monomial_order>)
        monomial_order::block_size = params.num_elim_vars + 1;

    monomial_base::exp_size =
        monomial_order::exponent_size(input_data.num_vars);
    monomial_base::initialize_weights(params.seed);

    basis_type basis{input_data.num_vars, input_data.field_char};
    basis.import_generators(input_data);

    GAMBA_DEBUG(basis.write_generators(debug);)

    if (params.basis_info)
        basis.print_info(std::cout);

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    f4_main(basis, params);

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.overall_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.overall_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    if (params.timings_info)
        print_timings();

    if (params.stats_info)
        print_statistics();

    generators_data output_data;
    basis.export_generators(output_data);

    /* copy variable names from input data */
    output_data.var_names = input_data.var_names;

    return output_data;
}

}  // namespace gamba
