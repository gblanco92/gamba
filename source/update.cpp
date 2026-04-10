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

#include "update.hpp"

namespace gamba
{

void update_f4(spair_set& spairs,
               base_polynomial_basis& basis,
               monomial_order const& mon_order,
               size_t const prev_num_gens)
{
    /* update process and specially G-M's speed is afected by virtual calls in
     * monomial_order abstract class, move back to templace metaprogramming */

    switch (mon_order.type())
    {
        case params::order::grevlex:
            update_f4_impl<order_grevlex>(spairs, basis, prev_num_gens);

            return;
        case params::order::deglex:
            update_f4_impl<order_deglex>(spairs, basis, prev_num_gens);

            return;
        case params::order::lexic:
            update_f4_impl<order_lexic>(spairs, basis, prev_num_gens);

            return;
        case params::order::blockelim:
            update_f4_impl<order_blockelim>(spairs, basis, prev_num_gens);

            return;
        case params::order::grevlexw:
            update_f4_impl<order_grevlexw>(spairs, basis, prev_num_gens);

            return;
    }
}

}  // namespace gamba
