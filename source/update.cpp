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

#include "order.hpp"

namespace gamba
{

#ifndef GAMBA_UNITY_BUILD
template void update_f4<order_grevlex>(spair_set& spairs,
                                       base_polynomial_basis& basis,
                                       size_t const prev_num_gens);
// template void update_f4<order_lexic>(spair_set& spairs,
//                                      base_polynomial_basis& basis,
//                                      size_t const prev_num_gens);
template void update_f4<order_blockelim>(spair_set& spairs,
                                         base_polynomial_basis& basis,
                                         size_t const prev_num_gens);
#endif

}  // namespace gamba
