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

#include <flint/ulong_extras.h>
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
    using integer_type = UIntT;

    // clang-format off

    template <class T> struct montgomery_{ using type = T; };
    template <> struct montgomery_<uint8_t> { using type = uint32_t; };
    template <> struct montgomery_<uint16_t> { using type = uint32_t; };
    template <> struct montgomery_<uint32_t> { using type = uint64_t; };

    // clang-format on

    using montgomery_type = montgomery_<integer_type>::type;

    /* if montgomery_type = uint32_t  --> r = 2^32;
       if montgomery_type = uint64_t  --> r = 2^64 */
    static constexpr size_t const r_nbits = CHAR_BIT * sizeof(montgomery_type);
    static constexpr auto const r_const = static_cast<uint128_t>(1) << r_nbits;

    /* computes the inverse of x mod 2^r_nbits */
    static constexpr montgomery_type inverse_mod2(montgomery_type const x)
    {
        constexpr size_t const r_log2_nbits = std::bit_width(r_nbits) - 1;

        montgomery_type xr = 1;
        for (size_t i = 0; i < r_log2_nbits; ++i)
            xr *= 2U - x * xr;

        return xr;
    }

    explicit field_traits(uint32_t const _n) :
            n{static_cast<integer_type>(_n)},

            r{static_cast<integer_type>(r_const % n)},

            r2{static_cast<integer_type>(n_powmod2(r, 2, n))},

            r4{static_cast<integer_type>(n_powmod2(r, 4, n))},

            /* do not use FLINT's ulong functions since FLINT_BITS < 63 */
            nr{inverse_mod2(n)}
    {
        assert(_n < r_const);
        assert((static_cast<uint64_t>(nr) * n) % r_const == 1);

        /* we add elements in the range [0, 2*n^2) while the sum is < n * r;
         * montgomery_max_fma >= 2^32 / 2 / 65521 = 32775 */
        size_t const montgomery_max_fma = r_const / 2U / n;

        size_t float_max_fma = std::numeric_limits<size_t>::max();  // infinity

        /* we add elements in the range [0, 2*n^2) while the sum is < 2^52;
         * float_max_fma >= 2^52 / 2 / 65521^2 = 524'528 */
        if constexpr (not std::is_same_v<integer_type, uint32_t>)
            float_max_fma = (1UL << 52) / (2U * n * n);
        /* AVX2 magic constant trick only works if < 2^52 */

        max_fma = std::min(montgomery_max_fma, float_max_fma);

        /* force low max_fma during testing */
        GAMBA_DEBUG(max_fma = 25);
    }

    integer_type add(integer_type const x, integer_type const y) const
    {
        auto const sum = n_addmod(x, y, n);

        return static_cast<integer_type>(sum);
    }

    integer_type multiply(integer_type const x, integer_type const y) const
    {
        auto const prod = n_mulmod2(x, y, n);

        return static_cast<integer_type>(prod);
    }

    integer_type inverse(integer_type const x) const
    {
        auto const inv = n_invmod(x, n);

        return static_cast<integer_type>(inv);
    }

    integer_type modular_reduce(mpz_class const& x) const
    {
        ulong const res = mpz_fdiv_ui(x.get_mpz_t(), n);

        return static_cast<integer_type>(res);
    }

    /* returns uint32_t to avoid overflow when 'n' has *exactly* 16 bits */
    uint32_t reduce(uint64_t const x) const
    {
        auto const res = montgomery<montgomery_type>::reduce(x, n, nr);

        /* WARNING: can overflow if coeff_type == uint32_t && n > 2^31 */
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

    /* r^4 modulo n */
    integer_type const r4;

    /* inverse of n modulo 2^r_num_bits */
    montgomery_type const nr;

    /* number of terms we can accumulate before a modular reduction */
    size_t max_fma;
};

template <>
struct field_traits<mpq_class>
{
    explicit field_traits([[maybe_unused]] uint32_t const n) { assert(n == 0); }

    static mpq_class add(mpq_class const& x, mpq_class const& y)
    {
        return x + y;
    }

    static mpq_class multiply(mpq_class const& x, mpq_class const& y)
    {
        return x * y;
    }

    static mpq_class inverse(mpq_class const& x) { return 1 / x; }
};

}  // namespace gamba
