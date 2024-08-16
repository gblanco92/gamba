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

#pragma once

#include <climits>

#include "utils.hpp"

namespace gamba
{

template <std::unsigned_integral IntT>
struct montgomery
{
    using integer_type  = IntT;
    using next_int_type = next_size_t<integer_type>;

    static constexpr size_t const num_bits = CHAR_BIT * sizeof(integer_type);

#ifdef __clang__
    __attribute__((no_sanitize("integer")))
#endif
    /* Prerequisite: x < r * n.
     * Returns a number in the [0, 2n) range */
    static integer_type
    reduce(next_int_type x, integer_type n, integer_type nr)
    {
        /* this overflows but its a feature not a bug */
        integer_type const q = static_cast<integer_type>(x) * nr;
        integer_type const m = static_cast<integer_type>(  // NOLINT
            (static_cast<next_int_type>(q) * n) >> num_bits);
        integer_type const y = static_cast<integer_type>(x >> num_bits) + n - m;

        return y;
    }

#ifdef __clang__
    __attribute__((no_sanitize("integer")))
#endif
    /* Prerequisite: x < r * n.
     * Returns a number in the [0, n) range */
    static integer_type
    reduce_normalize(next_int_type x, integer_type n, integer_type nr)
    {
        /* this overflows but its a feature not a bug */
        integer_type const q  = static_cast<integer_type>(x) * nr;
        next_int_type const m = static_cast<next_int_type>(q) * n;
        integer_type const y  = static_cast<integer_type>((x - m) >> num_bits);

        return x < m ? y + n : y;
    }
};

}  // namespace gamba
