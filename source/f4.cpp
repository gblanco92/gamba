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

#include "f4.hpp"

#include "logger.hpp"

namespace gamba
{

#ifndef GAMBA_UNITY_BUILD
extern template class polynomial_basis<uint8_t>;
extern template class polynomial_basis<uint16_t>;
extern template class polynomial_basis<uint32_t>;

extern template class linalg_v3<uint8_t>;
extern template class linalg_v3<uint16_t>;
extern template class linalg_v3<uint32_t>;

extern template void update_f4<order_grevlex>(spair_set& spairs,
                                              base_polynomial_basis& basis,
                                              size_t const prev_num_gens);
// extern template void update_f4<order_lexic>(spair_set& spairs,
//                                             base_polynomial_basis& basis,
//                                             size_t const prev_num_gens);
extern template void update_f4<order_blockelim>(spair_set& spairs,
                                                base_polynomial_basis& basis,
                                                size_t const prev_num_gens);

extern template class std::pair<size_t, size_t> select_spairs<order_grevlex>(
    spair_set& spairs);
// extern template class std::pair<size_t, size_t> select_spairs<order_lexic>(
//     spair_set& spairs);
extern template class std::pair<size_t, size_t> select_spairs<order_blockelim>(
    spair_set& spairs);

extern template void reduce<order_grevlex, uint8_t>(
    polynomial_basis<uint8_t>& basis,
    matrix_f4& matrix);
extern template void reduce<order_grevlex, uint16_t>(
    polynomial_basis<uint16_t>& basis,
    matrix_f4& matrix);
extern template void reduce<order_grevlex, uint32_t>(
    polynomial_basis<uint32_t>& basis,
    matrix_f4& matrix);

// extern template void reduce<order_lexic, uint8_t>(
//      polynomial_basis<uint8_t>& basis,
//      matrix_f4& matrix);
// extern template void reduce<order_lexic, uint16_t>(
//      polynomial_basis<uint16_t>& basis,
//      matrix_f4& matrix);
// extern template void reduce<order_lexic, uint32_t>(
//      polynomial_basis<uint32_t>& basis,
//      matrix_f4& matrix);

extern template void reduce<order_blockelim, uint8_t>(
    polynomial_basis<uint8_t>& basis,
    matrix_f4& matrix);
extern template void reduce<order_blockelim, uint16_t>(
    polynomial_basis<uint16_t>& basis,
    matrix_f4& matrix);
extern template void reduce<order_blockelim, uint32_t>(
    polynomial_basis<uint32_t>& basis,
    matrix_f4& matrix);

template void f4_main<order_grevlex, uint8_t>(
    polynomial_basis<uint8_t>& basis);  // NOFORMAT
template void f4_main<order_grevlex, uint16_t>(
    polynomial_basis<uint16_t>& basis);
template void f4_main<order_grevlex, uint32_t>(
    polynomial_basis<uint32_t>& basis);

// template void f4_main<order_lexic, uint8_t>(
//   polynomial_basis<uint8_t>& basis);
// template void f4_main<order_lexic, uint16_t>(
//   polynomial_basis<uint16_t>& basis);
// template void f4_main<order_lexic, uint32_t>(
//   polynomial_basis<uint32_t>& basis);

template void f4_main<order_blockelim, uint8_t>(
    polynomial_basis<uint8_t>& basis);
template void f4_main<order_blockelim, uint16_t>(
    polynomial_basis<uint16_t>& basis);
template void f4_main<order_blockelim, uint32_t>(
    polynomial_basis<uint32_t>& basis);
#endif

namespace f4
{

void print_column_names()
{
    log::print(log::INFO2, "\n┌{0:─^118}┐\n", "");

    log::print(log::INFO2, "│ {}{:>12}{:>19}{:>15}{:>20}{:>19}{:>14}{:>17}",
               "deg", "spairs", "matrix", "density", "new gens", "ech. time",
               "mem. usage", "total time │\n");

    log::print(log::INFO2, "├{:─^118}┤\n", "");
}

void print_bottom_line()
{
    if (params::no_reduce)
        log::print(log::INFO2, "╘{:═^118}╛\n", "");
    else
        log::print(log::INFO2, "╞{:═^118}╡\n", "");
}

}  // namespace f4

}  // namespace gamba
