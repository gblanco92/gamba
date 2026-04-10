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
extern template void f4_main<uint8_t>(polynomial_basis<uint8_t>& basis,
                                      monomial_order const& mon_order,
                                      learn_f4_data* learn_data);
extern template void f4_main<uint16_t>(polynomial_basis<uint16_t>& basis,
                                       monomial_order const& mon_order,
                                       learn_f4_data* learn_data);
extern template void f4_main<uint32_t>(polynomial_basis<uint32_t>& basis,
                                       monomial_order const& mon_order,
                                       learn_f4_data* learn_data);

template generators_data groebner_basis_main<uint8_t>(
    generators_data const& input_data,
    monomial_order const& mon_order);
template generators_data groebner_basis_main<uint16_t>(
    generators_data const& input_data,
    monomial_order const& mon_order);
template generators_data groebner_basis_main<uint32_t>(
    generators_data const& input_data,
    monomial_order const& mon_order);
template generators_data groebner_basis_main<fmpq_class>(
    generators_data const& input_data,
    monomial_order const& mon_order);
#endif

}  // namespace gamba
