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

extern template void reduce<uint8_t>(polynomial_basis<uint8_t>& basis,
                                     monomial_order const& mon_basis,
                                     matrix_f4& matrix);
extern template void reduce<uint16_t>(polynomial_basis<uint16_t>& basis,
                                      monomial_order const& mon_basis,
                                      matrix_f4& matrix);
extern template void reduce<uint32_t>(polynomial_basis<uint32_t>& basis,
                                      monomial_order const& mon_basis,
                                      matrix_f4& matrix);

template void f4_main<uint8_t>(polynomial_basis<uint8_t>& basis,
                               monomial_order const& mon_order,
                               learn_f4_data* learn_data);
template void f4_main<uint16_t>(polynomial_basis<uint16_t>& basis,
                                monomial_order const& mon_order,
                                learn_f4_data* learn_data);
template void f4_main<uint32_t>(polynomial_basis<uint32_t>& basis,
                                monomial_order const& mon_order,
                                learn_f4_data* learn_data);
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

void print_time(std::chrono::duration<double> const time, bool const new_line)
{
    log::print(log::INFO2, "{:10.2f} sec", time.count());

    if (new_line)
        log::print(log::INFO2, " │\n");

    ::fflush(stdout);
}

}  // namespace f4

}  // namespace gamba
