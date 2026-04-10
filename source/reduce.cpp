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

#include "reduce.hpp"

namespace gamba
{

#ifndef GAMBA_UNITY_BUILD
extern template class polynomial_basis<uint8_t>;
extern template class polynomial_basis<uint16_t>;
extern template class polynomial_basis<uint32_t>;

extern template class linalg<uint8_t>;
extern template class linalg<uint16_t>;
extern template class linalg<uint32_t>;

template void reduce<uint8_t>(polynomial_basis<uint8_t>& basis,
                              monomial_order const& mon_order,
                              matrix_f4& matrix);
template void reduce<uint16_t>(polynomial_basis<uint16_t>& basis,
                               monomial_order const& mon_order,
                               matrix_f4& matrix);
template void reduce<uint32_t>(polynomial_basis<uint32_t>& basis,
                               monomial_order const& mon_order,
                               matrix_f4& matrix);
#endif

}  // namespace gamba
