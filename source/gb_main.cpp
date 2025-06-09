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

#include "gb_main.hpp"

namespace gamba
{

#ifndef GAMBA_UNITY_BUILD
extern template void f4_main<order_grevlex, uint8_t>(
    polynomial_basis<uint8_t>& basis);
extern template void f4_main<order_grevlex, uint16_t>(
    polynomial_basis<uint16_t>& basis);
extern template void f4_main<order_grevlex, uint32_t>(
    polynomial_basis<uint32_t>& basis);

// extern template void f4_main<order_lexic, uint8_t>(
//   polynomial_basis<uint8_t>& basis);
// extern template void f4_main<order_lexic, uint16_t>(
//   polynomial_basis<uint16_t>& basis);
// extern template void f4_main<order_lexic, uint32_t>(
//   polynomial_basis<uint32_t>& basis);

extern template void f4_main<order_blockelim, uint8_t>(
    polynomial_basis<uint8_t>& basis);
extern template void f4_main<order_blockelim, uint16_t>(
    polynomial_basis<uint16_t>& basis);
extern template void f4_main<order_blockelim, uint32_t>(
    polynomial_basis<uint32_t>& basis);

extern template void f4_modular<order_grevlex>(
    polynomial_basis<mpq_class>& basis);
// extern template void f4_modular<order_lexic>(
//     polynomial_basis<mpq_class>& basis);
extern template void f4_modular<order_blockelim>(
    polynomial_basis<mpq_class>& basis);

template generators_data groebner_basis_main<uint8_t, order_grevlex>(
    generators_data const& input_data);
template generators_data groebner_basis_main<uint16_t, order_grevlex>(
    generators_data const& input_data);
template generators_data groebner_basis_main<uint32_t, order_grevlex>(
    generators_data const& input_data);
template generators_data groebner_basis_main<mpq_class, order_grevlex>(
    generators_data const& input_data);

// template generators_data groebner_basis_main<uint8_t, order_lexic>(
//  generators_data const& input_data);
// template generators_data groebner_basis_main<uint16_t, order_lexic>(
//  generators_data const& input_data);
// template generators_data groebner_basis_main<uint32_t, order_lexic>(
//  generators_data const& input_data);
// template generators_data groebner_basis_main<mpq_class, order_lexic>(
//  generators_data const& input_data);

template generators_data groebner_basis_main<uint8_t, order_blockelim>(
    generators_data const& input_data);
template generators_data groebner_basis_main<uint16_t, order_blockelim>(
    generators_data const& input_data);
template generators_data groebner_basis_main<uint32_t, order_blockelim>(
    generators_data const& input_data);
template generators_data groebner_basis_main<mpq_class, order_blockelim>(
    generators_data const& input_data);
#endif

}  // namespace gamba
