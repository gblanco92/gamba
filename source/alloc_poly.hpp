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

#pragma once

#include <utility>

#include "monomial.hpp"

namespace gamba
{

using basis_monomial  = monomial<basis_hashtable>;
using spair_monomial  = monomial<spair_hashtable>;
using matrix_monomial = monomial<matrix_hashtable>;

template <class CoefficientType>
std::pair<basis_monomial*, CoefficientType*> allocate_polynomial(
    size_t const num_terms)
{
    using coeff_type    = CoefficientType;
    using monomial_type = basis_monomial;

    constexpr size_t const std_align = __STDCPP_DEFAULT_NEW_ALIGNMENT__;

    size_t const num_bytes_mons  = num_terms * sizeof(monomial_type);
    size_t const num_bytes_coefs = num_terms * sizeof(coeff_type);
    size_t const num_bytes_padd  = padding<std_align>(num_bytes_coefs);

    size_t const num_bytes = num_bytes_mons + num_bytes_coefs + num_bytes_padd;
    size_t const offset_bytes = num_bytes_coefs + num_bytes_padd;

    /* allocate monomials & coefficients in a single allocation; deallocate
     * memory by deallocating the coefficients; use ::operator new/delete to
     * avoid casting coefficients pointers back to char when deallocating */
    void* const mem_ptr = ::operator new(num_bytes);
    auto* const cfs_ptr = reinterpret_cast<coeff_type*>(mem_ptr);
    auto* const mon_ptr = reinterpret_cast<monomial_type*>(mem_ptr)
                        + offset_bytes / sizeof(monomial_type);

    /* both memory blocks have the same alignment as if allocated separately */
    assert(reinterpret_cast<std::uintptr_t>(mem_ptr) % std_align == 0);
    assert(reinterpret_cast<std::uintptr_t>(cfs_ptr) % std_align == 0);

    return std::make_pair(mon_ptr, cfs_ptr);
}

}  // namespace gamba
