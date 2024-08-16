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

#include <cassert>
#include <climits>
#include <concepts>
#include <type_traits>

#include <gmpxx.h>

#include "montgomery.hpp"
#include "utils.hpp"

namespace gamba
{

template <class CoefficientType>
struct field_traits;

template <std::unsigned_integral UIntT>
struct field_traits<UIntT>
{
    using integer_type  = UIntT;
    using next_int_type = next_size_t<integer_type>;

    /* montgomery_type is either uint32_t or uint64_t */
    using montgomery_type =
        std::conditional_t<std::is_same_v<integer_type, uint32_t>,
                           uint64_t,
                           uint32_t>;

    /* if integer_type = uint16_t or uint8_t --> r = 2^32;
       if integer_type = uint32_t            --> r = 2^64 */
    static constexpr size_t const r_nbits = CHAR_BIT * sizeof(montgomery_type);
    static constexpr size_t const r_log2_nbits = std::bit_width(r_nbits) - 1;
    static constexpr auto const r_const = static_cast<uint128_t>(1) << r_nbits;

    /* computes the inverse of x mod 2^r_nbits */
    static constexpr montgomery_type inverse_mod2(montgomery_type const x)
    {
        montgomery_type xr{1};
        for (size_t i = 0; i < r_log2_nbits; ++i)
        {
            xr *= 2U - x * xr;
        }

        return xr;
    }

    explicit field_traits(uint32_t const _n) :
            n{static_cast<integer_type>(_n)},

            r{static_cast<integer_type>(r_const % n)},

            r2{static_cast<integer_type>(static_cast<next_int_type>(r) * r
                                         % n)},

            r3{static_cast<integer_type>(static_cast<next_int_type>(r2) * r
                                         % n)},

            r4{static_cast<integer_type>(static_cast<next_int_type>(r2) * r2
                                         % n)},

            nr{inverse_mod2(n)},

            /* we are adding elements in the range [0, 2*n^2) while the sum is
               < n * r, hence the sum can have at most (r / 2 / n) terms */
            max_fma{static_cast<size_t>(r_const / 2U / n)}
    {
        assert(_n < r_const);
        assert((static_cast<uint64_t>(nr) * n) % r_const == 1);
    }

    /* computes the inverse of x mod n */
    integer_type inverse(integer_type const _x) const
    {
        using signed_next_type = std::make_signed_t<next_int_type>;

        assert(_x % n != 0);

        signed_next_type t;
        signed_next_type q;
        signed_next_type x0 = 0;
        signed_next_type x1 = 1;

        auto n1 = static_cast<signed_next_type>(n);
        auto x  = static_cast<signed_next_type>(_x);

        while (x > 1)
        {
            q  = x / n1;
            t  = n1;
            n1 = x % n1;
            x  = t;
            t  = x0;
            x0 = x1 - static_cast<signed_next_type>(q * x0);
            x1 = t;
        }

        if (x1 < 0)
        {
            x1 += static_cast<signed_next_type>(n);
        }

        return static_cast<integer_type>(x1);
    }

    integer_type add(integer_type const x, integer_type const y) const
    {
        auto const sum =
            (static_cast<next_int_type>(x) + y) % static_cast<next_int_type>(n);

        return static_cast<integer_type>(sum);
    }

    integer_type multiply(integer_type const x, integer_type const y) const
    {
        auto const prod =
            (static_cast<next_int_type>(x) * y) % static_cast<next_int_type>(n);

        return static_cast<integer_type>(prod);
    }

    /* Return uint32_t to avoid overflow when n has exactly 16 bits */
    uint32_t reduce(uint64_t const x) const
    {
        auto const res = montgomery<montgomery_type>::reduce(x, n, nr);

        /* WARNING: this can overflow if coeff_type == uint32_t && n > 2^31 */
        return static_cast<uint32_t>(res);
    }

    integer_type reduce_normalize(uint64_t const x) const
    {
        auto const res =
            montgomery<montgomery_type>::reduce_normalize(x, n, nr);

        return static_cast<integer_type>(res);
    }

    integer_type transform(uint64_t const x) const
    {
        auto const res =
            montgomery<montgomery_type>::reduce_normalize(x * r2, n, nr);

        return static_cast<integer_type>(res);
    }

    /* field characteristic */
    integer_type const n;

    /* 2^r_num_bits modulo n */
    integer_type const r;

    /* r^2 modulo n */
    integer_type const r2;

    /* r^3 modulo n */
    integer_type const r3;

    /* r^4 modulo n */
    integer_type const r4;

    /* inverse of n modulo 2^r_num_bits */
    montgomery_type const nr;

    /* number of terms we can accumulate before a modular reduction */
    size_t const max_fma;
};

template <>
struct field_traits<mpz_class>
{
    explicit field_traits([[maybe_unused]] uint32_t const n) { assert(n == 0); }

    static mpz_class add(mpz_class const& x, mpz_class const& y)
    {
        return x + y;
    }

    static mpz_class multiply(mpz_class const& x, mpz_class const& y)
    {
        return x * y;
    }
};

}  // namespace gamba
