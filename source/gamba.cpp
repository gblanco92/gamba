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

#include <stdexcept>

#include "gamba.hpp"

namespace gamba
{

generators_data groebner_basis(generators_data const& data,
                               gamba_params& params)
{
    switch (params.mon_order)
    {
        case order_options::grevlex:
            return groebner_basis_coeff_type<order_grevlex>(data, params);
#if 0
        case order_options::lexic:
            return groebner_basis_coeff_type<order_lexic>(data, params);
#endif
        case order_options::blockelim:
            return groebner_basis_coeff_type<order_blockelim>(data, params);
        default:
            throw std::runtime_error{"Monomial order not supported."};
    }
}

}  // namespace gamba
