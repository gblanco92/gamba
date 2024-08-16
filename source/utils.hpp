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

#include <algorithm>
#include <concepts>
#include <numeric>
#include <vector>

#include "config.hpp"

#ifdef DEBUG
#    define FORCE_INLINE inline
#else
#    define FORCE_INLINE inline __attribute__((__always_inline__))
#endif

#define LIKELY(x)   __builtin_expect_with_probability((x), 1, 1.0)
#define UNLIKELY(x) __builtin_expect_with_probability((x), 0, 1.0)

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(X)        STRINGIFY_HELPER(X)

#ifdef __clang__
#    define PRAGMA_UNROLL_LOOP(x) _Pragma(STRINGIFY(clang loop unroll_count(x)))
#else
#    define PRAGMA_UNROLL_LOOP(x) _Pragma(STRINGIFY(GCC unroll(x)))
#endif

namespace gamba
{

/* some random STL-like algorithm functions */

template <class T, class Compare>
[[nodiscard]] std::vector<size_t> sort_permutation(std::vector<T> const& v,
                                                   Compare comp)
{
    std::vector<size_t> p(v.size());
    std::iota(std::begin(p), std::end(p), 0);
    std::sort(std::begin(p), std::end(p),
              [&](size_t const i, size_t const j) { return comp(v[i], v[j]); });

    return p;
}

template <class T>
[[nodiscard]] std::vector<T> apply_permutation(std::vector<T> const& v,
                                               std::vector<size_t> const& p)
{
    std::vector<T> s(v.size());
    std::transform(std::cbegin(p), std::cend(p), std::begin(s),
                   [&](size_t const i) { return v[i]; });

    return s;
}

/* helper functions to operate safely with mod p coefficients */

using int128_t  = __int128_t;
using uint128_t = __uint128_t;

template <class>
struct next_size;
template <class T>
using next_size_t = typename next_size<T>::type;

template <class T>
struct tag
{
    using type = T;
};

template <>
struct next_size<int8_t> : tag<int16_t>
{};

template <>
struct next_size<int16_t> : tag<int32_t>
{};

template <>
struct next_size<int32_t> : tag<int64_t>
{};

template <>
struct next_size<int64_t> : tag<int128_t>
{};

template <>
struct next_size<uint8_t> : tag<uint16_t>
{};

template <>
struct next_size<uint16_t> : tag<uint32_t>
{};

template <>
struct next_size<uint32_t> : tag<uint64_t>
{};

template <>
struct next_size<uint64_t> : tag<uint128_t>
{};

/* helper definitions for accessing tuple elements with better namings */

template <size_t N>
struct get_t
{
    template <class T>
    decltype(auto) operator()(T&& t) const
    {
        return std::get<N>(std::forward<T>(t));
    }
};

inline constexpr get_t<0> get_mon;
inline constexpr get_t<1> get_coef;

inline constexpr get_t<0> get_mons;
inline constexpr get_t<1> get_coefs;
inline constexpr get_t<2> get_deg;

/* round up the number 'num' to the closest multiple of 'mult' */
template <std::unsigned_integral UIntT>
constexpr UIntT round_up(UIntT num, UIntT mult)
{
    assert(mult);
    return ((num + mult - 1) / mult) * mult;
}

/* round down the number 'num' to the closest multiple of 'mult' */
template <std::unsigned_integral UIntT>
constexpr UIntT round_down(UIntT num, UIntT mult)
{
    assert(mult);
    return (num / mult) * mult;
}

/* number of padding bytes needed for 'size' to be a multiple of Mult */
template <size_t Mult>
constexpr size_t padding(size_t size)
{
    return ((Mult - size) % Mult) % Mult;
}

/* non-owning wrapper for a buffer row (strided memory)*/
template <class ValueType, class IndexType, size_t StrideOffset>
struct buffer_row_wrapper
{
    using value_type = ValueType;
    using index_type = IndexType;

    bool operator>(buffer_row_wrapper const& other) const
    {
        return *this->m_idx > *other.m_idx;
    }

    static inline size_t size;

    value_type* m_ptr;
    index_type* m_idx;
};

/* swap two (strided) rows of a linear algebra buffer in place */
template <class ValueType>
void swap_buffer_rows(ValueType* lhs_buffer,
                      size_t const lhs_idx,
                      ValueType* rhs_buffer,
                      size_t const rhs_idx,
                      size_t const buffer_size)
{
    size_t i = 0;

    for (; i < (buffer_size & ~0x3U); i += 4)
    {
        std::swap(lhs_buffer[lhs_idx + 0 * NUM_ROWS_BUFFER],
                  rhs_buffer[rhs_idx + 0 * NUM_ROWS_BUFFER]);

        std::swap(lhs_buffer[lhs_idx + 1 * NUM_ROWS_BUFFER],
                  rhs_buffer[rhs_idx + 1 * NUM_ROWS_BUFFER]);

        std::swap(lhs_buffer[lhs_idx + 2 * NUM_ROWS_BUFFER],
                  rhs_buffer[rhs_idx + 2 * NUM_ROWS_BUFFER]);

        std::swap(lhs_buffer[lhs_idx + 3 * NUM_ROWS_BUFFER],
                  rhs_buffer[rhs_idx + 3 * NUM_ROWS_BUFFER]);

        lhs_buffer += 4 * NUM_ROWS_BUFFER;
        rhs_buffer += 4 * NUM_ROWS_BUFFER;
    }

    for (; i < buffer_size; ++i)
    {
        std::swap(lhs_buffer[lhs_idx], rhs_buffer[rhs_idx]);

        lhs_buffer += NUM_ROWS_BUFFER;
        rhs_buffer += NUM_ROWS_BUFFER;
    }
}

template <class ValueType, class IndexType, size_t StrideOffset>
void swap(buffer_row_wrapper<ValueType, IndexType, StrideOffset>& lhs,
          buffer_row_wrapper<ValueType, IndexType, StrideOffset>& rhs) noexcept
{
    using row_type = buffer_row_wrapper<ValueType, IndexType, StrideOffset>;

    swap_buffer_rows(lhs.m_ptr, 0, rhs.m_ptr, 0, row_type::size);

    std::iter_swap(lhs.m_idx, rhs.m_idx);
}

template <class T, class Allocator>
double memory_size(std::vector<T, Allocator> const& v)
{
    /* return memory footprint of vector in megabytes */
    return static_cast<double>(v.capacity() * sizeof(T)) / 1024.0 / 1024.0;
}

template <class T>
double memory_size(std::span<T> const& v)
{
    /* return memory footprint of span in megabytes */
    return static_cast<double>(v.size() * sizeof(T)) / 1024.0 / 1024.0;
}

}  // namespace gamba
