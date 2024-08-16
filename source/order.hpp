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

#include <type_traits>

namespace gamba
{

/* degree reverse lexicographical order */
struct order_grevlex
{
    using degree_order  = void;
    using reverse_order = void;

    /* comparators return a signed type so negation for reverse sorting works */
    template <class MonomialType>
    constexpr int32_t operator()(MonomialType const& lhs,
                                 MonomialType const& rhs) const
    {
        auto const degl = degree(lhs);
        auto const degr = degree(rhs);

        if (degl < degr)
            return -1;
        if (degl != degr)
            return 1;

        auto const* const expl = lhs.cbegin();
        auto const* const expr = rhs.cbegin();

        size_t i;
        for (i = lhs.size() - 1; i > 1 and expr[i] == expl[i]; --i)
            ;

        return static_cast<int32_t>(expr[i]) - static_cast<int32_t>(expl[i]);
    }

    template <class MonomialType>
    constexpr static auto degree(MonomialType const& mon)
        -> MonomialType::exponent_type
    {
        return mon.cbegin()[0];
    }

    constexpr static size_t exponent_size(size_t const num_vars)
    {
        return num_vars + 1;
    }

    static inline size_t block_size{0};
};

/* lexicographical order */
struct order_lexic
{
    template <class MonomialType>
    constexpr int32_t operator()(MonomialType const& lhs,
                                 MonomialType const& rhs) const
    {
        auto const* const expl = lhs.cbegin();
        auto const* const expr = rhs.cbegin();

        size_t i;
        for (i = 1; i < lhs.size() - 1 and expl[i] == expr[i]; ++i)
            ;

        return static_cast<int32_t>(expl[i]) - static_cast<int32_t>(expr[i]);
    }

    template <class MonomialType>
    constexpr static auto degree(MonomialType const& mon)
        -> MonomialType::exponent_type
    {
        return mon.cbegin()[0];
    }

    constexpr static size_t exponent_size(size_t const num_vars)
    {
        return num_vars + 1;
    }

    static inline size_t block_size{0};
};

/* two blocks degree reverse lexicographical order */
struct order_blockelim
{
    using degree_order  = void;
    using reverse_order = void;
    using block_order   = void;

    template <class MonomialType>
    constexpr int32_t operator()(MonomialType const& lhs,
                                 MonomialType const& rhs) const
    {
        auto const* const expl = lhs.cbegin();
        auto const* const expr = rhs.cbegin();

        /* first block */
        auto const deg1l = expl[0];
        auto const deg1r = expr[0];

        if (deg1l < deg1r)
            return -1;
        if (deg1l != deg1r)
            return 1;

        size_t i;
        for (i = block_size - 1; i > 0 and expl[i] == expr[i]; --i)
            ;

        if (i != 0)
        {
            return static_cast<int32_t>(expr[i])
                 - static_cast<int32_t>(expl[i]);
        }

        /* second block */
        auto const deg2l = expl[block_size];
        auto const deg2r = expr[block_size];

        if (deg2l < deg2r)
            return -1;
        if (deg2l != deg2r)
            return 1;

        for (i = lhs.size() - 1; i > 1 and expl[i] == expr[i]; --i)
            ;

        return static_cast<int32_t>(expr[i]) - static_cast<int32_t>(expl[i]);
    }

    template <class MonomialType>
    constexpr static auto degree(MonomialType const& mon)
        -> MonomialType::exponent_type
    {
        auto const* const exp = mon.cbegin();
        return (exp[0] + exp[block_size]);
    }

    constexpr static size_t exponent_size(size_t const num_vars)
    {
        return num_vars + 2;
    }

    /* number of elimination variables + degree of the block */
    static inline size_t block_size{};
};

template <class, class = void>
struct is_degree_order : std::false_type
{};

/* specialization recognizes monomial order that have nested ::degree_order */
template <class T>
struct is_degree_order<T, std::void_t<typename T::degree_order>>
        : std::true_type
{};

template <class T>
inline constexpr bool is_degree_order_v = is_degree_order<T>::value;

template <class, class = void>
struct is_reverse_order : std::false_type
{};

/* specialization recognizes monomial order that have nested ::reverse_order */
template <class T>
struct is_reverse_order<T, std::void_t<typename T::reverse_order>>
        : std::true_type
{};

template <class T>
inline constexpr bool is_reverse_order_v = is_reverse_order<T>::value;

template <class, class = void>
struct is_block_order : std::false_type
{};

/* specialization recognizes monomial order that have nested ::block_order */
template <class T>
struct is_block_order<T, std::void_t<typename T::block_order>> : std::true_type
{};

template <class T>
inline constexpr bool is_block_order_v = is_block_order<T>::value;

}  // namespace gamba
