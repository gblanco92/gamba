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

#include <limits>
#include <stdexcept>

#include "gamba.hpp"

#include "logger.hpp"

namespace gamba
{

template <class CoefficientType, class MonomialOrder>
generators_data groebner_basis_main(generators_data const& data);  // NOLINT

template <class CoefficientType>
generators_data groebner_basis_mon_order(generators_data const& data)  // NOLINT
{
    using coeff_type = CoefficientType;

    switch (params::mon_order)
    {
        case params::order::grevlex:
            return groebner_basis_main<coeff_type, order_grevlex>(data);
#if 0
        case params::order::lexic:
            return groebner_basis_main<coeff_type, order_lexic>(data, params);
#endif
        case params::order::blockelim:
            return groebner_basis_main<coeff_type, order_blockelim>(data);
        default:
            throw std::runtime_error{"Monomial order not supported."};
    }
}

generators_data groebner_basis(generators_data const& data)
{
    if (data.field_char == 2)
    {
        throw std::runtime_error("Field characteristic = 2 not supported.");
    }

    if (data.field_char == 0)
    {
        throw std::runtime_error{"Rational numbers not supported yet."};

        return groebner_basis_mon_order<mpq_class>(data);
    }

    if (data.field_char < std::numeric_limits<uint8_t>::max())
    {
        try
        {
            return groebner_basis_mon_order<uint8_t>(data);
        }
        catch (std::overflow_error const&)
        {
            log::print(log::WARN, "Index overflow detected. Restarting...\n");

            return groebner_basis_mon_order<uint32_t>(data);
        }
    }

    if (data.field_char < std::numeric_limits<uint16_t>::max())
    {
        try
        {
            return groebner_basis_mon_order<uint16_t>(data);
        }
        catch (std::overflow_error const&)
        {
            log::print(log::WARN, "Index overflow detected. Restarting...\n");

            return groebner_basis_mon_order<uint32_t>(data);
        }
    }

    /* max. allowed field characteristic is 2^31 - 1 */
    if (data.field_char <= std::numeric_limits<int32_t>::max())
    {
        return groebner_basis_mon_order<uint32_t>(data);
    }

    throw std::runtime_error("Field characteristic too large.");
}

}  // namespace gamba
