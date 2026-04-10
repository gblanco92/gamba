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

#include "monomial.hpp"
#include "params.hpp"

namespace gamba
{

/* virtual base monomial_order class */
struct monomial_order
{
    using degree_type = uint32_t;

    virtual ~monomial_order() = default;

    virtual int32_t cmp(monomial<basis_hashtable> const lhs,
                        monomial<basis_hashtable> const rhs) const = 0;

    virtual int32_t cmp(monomial<spair_hashtable> const lhs,
                        monomial<spair_hashtable> const rhs) const = 0;

    virtual int32_t cmp(monomial<matrix_hashtable> const lhs,
                        monomial<matrix_hashtable> const rhs) const = 0;

    virtual degree_type deg(monomial<basis_hashtable> const mon) const = 0;

    virtual degree_type deg(monomial<spair_hashtable> const mon) const = 0;

    virtual degree_type deg(monomial<matrix_hashtable> const mon) const = 0;

    virtual size_t exponent_size(size_t const num_vars) const = 0;

    virtual bool is_degree_order() const   = 0;
    virtual bool is_reverse_order() const  = 0;
    virtual bool is_block_order() const    = 0;
    virtual bool is_weighted_order() const = 0;

    virtual params::order type() const = 0;
};

template <class T>
struct monomial_order_interface : public monomial_order
{
private:
    monomial_order_interface() = default;

public:
    int32_t cmp(monomial<basis_hashtable> const lhs,
                monomial<basis_hashtable> const rhs) const final
    {
        return static_cast<T const*>(this)->operator()(lhs, rhs);
    }

    int32_t cmp(monomial<spair_hashtable> const lhs,
                monomial<spair_hashtable> const rhs) const final
    {
        return static_cast<T const*>(this)->operator()(lhs, rhs);
    }

    int32_t cmp(monomial<matrix_hashtable> const lhs,
                monomial<matrix_hashtable> const rhs) const final
    {
        return static_cast<T const*>(this)->operator()(lhs, rhs);
    }

    template <class MonomialType>
    constexpr static degree_type degree(MonomialType const mon)
    {
        return mon.cbegin()[0];
    }

    constexpr size_t exponent_size(size_t const num_vars) const override
    {
        return num_vars + 1;
    }

    degree_type deg(monomial<basis_hashtable> const mon) const final
    {
        return static_cast<T const*>(this)->degree(mon);
    }

    degree_type deg(monomial<spair_hashtable> const mon) const final
    {
        return static_cast<T const*>(this)->degree(mon);
    }

    degree_type deg(monomial<matrix_hashtable> const mon) const final
    {
        return static_cast<T const*>(this)->degree(mon);
    }

    bool is_degree_order() const override { return false; }

    bool is_reverse_order() const override { return false; }

    bool is_block_order() const override { return false; }

    bool is_weighted_order() const override { return false; }

    friend T;
};

/* degree reverse lexicographical order */
struct order_grevlex : monomial_order_interface<order_grevlex>
{
    using degree_order  = void;
    using reverse_order = void;

    /* comparators return a signed type so negation for reverse sorting works */
    template <class MonomialType>
    constexpr static int32_t operator()(MonomialType const lhs,
                                        MonomialType const rhs)
    {
        degree_type const degl = degree(lhs);
        degree_type const degr = degree(rhs);

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

    bool is_degree_order() const final { return true; }

    bool is_reverse_order() const final { return true; }

    params::order type() const override { return params::order::grevlex; }
};

/* degree lexicographical order */
struct order_deglex : monomial_order_interface<order_deglex>
{
    using degree_order = void;

    /* comparators return a signed type so negation for reverse sorting works */
    template <class MonomialType>
    constexpr static int32_t operator()(MonomialType const lhs,
                                        MonomialType const rhs)
    {
        degree_type const degl = degree(lhs);
        degree_type const degr = degree(rhs);

        if (degl < degr)
            return -1;
        if (degl != degr)
            return 1;

        auto const* const expl = lhs.cbegin();
        auto const* const expr = rhs.cbegin();

        size_t i;
        for (i = 1; i < lhs.size() - 1 and expl[i] == expr[i]; ++i)
            ;

        return static_cast<int32_t>(expl[i]) - static_cast<int32_t>(expr[i]);
    }

    bool is_degree_order() const final { return true; }

    params::order type() const final { return params::order::deglex; }
};

/* lexicographical order */
struct order_lexic : public monomial_order_interface<order_lexic>
{
    template <class MonomialType>
    constexpr static int32_t operator()(MonomialType const lhs,
                                        MonomialType const rhs)
    {
        auto const* const expl = lhs.cbegin();
        auto const* const expr = rhs.cbegin();

        size_t i;
        for (i = 1; i < lhs.size() - 1 and expl[i] == expr[i]; ++i)
            ;

        return static_cast<int32_t>(expl[i]) - static_cast<int32_t>(expr[i]);
    }

    params::order type() const final { return params::order::lexic; }
};

/* two blocks degree reverse lexicographical order */
struct order_blockelim : public monomial_order_interface<order_blockelim>
{
    using degree_order  = void;
    using reverse_order = void;
    using block_order   = void;

    template <class MonomialType>
    constexpr static int32_t operator()(MonomialType const lhs,
                                        MonomialType const rhs)
    {
        using monomial_type = MonomialType;

        size_t const ebz = monomial_type::block_size;

        auto const* const expl = lhs.cbegin();
        auto const* const expr = rhs.cbegin();

        /* first block */
        degree_type const deg1l = expl[0];
        degree_type const deg1r = expr[0];

        if (deg1l < deg1r)
            return -1;
        if (deg1l != deg1r)
            return 1;

        size_t i;
        for (i = ebz - 1; i > 0 and expl[i] == expr[i]; --i)
            ;

        if (i != 0)
        {
            return static_cast<int32_t>(expr[i])
                 - static_cast<int32_t>(expl[i]);
        }

        /* second block */
        auto const deg2l = expl[ebz];
        auto const deg2r = expr[ebz];

        if (deg2l < deg2r)
            return -1;
        if (deg2l != deg2r)
            return 1;

        for (i = lhs.size() - 1; i > 1 and expl[i] == expr[i]; --i)
            ;

        return static_cast<int32_t>(expr[i]) - static_cast<int32_t>(expl[i]);
    }

    template <class MonomialType>
    constexpr static degree_type degree(MonomialType const mon)
    {
        auto const* const exp = mon.cbegin();
        return (exp[0] + exp[MonomialType::block_size]);
    }

    constexpr size_t exponent_size(size_t const num_vars) const override
    {
        return num_vars + 2;
    }

    bool is_degree_order() const final { return true; }

    bool is_reverse_order() const final { return true; }

    bool is_block_order() const final { return true; }

    params::order type() const final { return params::order::blockelim; }
};

/* degree reverse weighted lexicographical order */
struct order_grevlexw : monomial_order_interface<order_grevlexw>
{
    using degree_order   = void;
    using reverse_order  = void;
    using weighted_order = void;

    template <class MonomialType>
    constexpr static degree_type degree(MonomialType const mon)
    {
        auto const* const exp = mon.cbegin();

        degree_type deg{0};
        for (size_t i = 1; i < mon.size(); ++i)
            deg += w_[i - 1] * exp[i];

        return deg;
    }

    template <class MonomialType>
    constexpr static int32_t operator()(MonomialType const lhs,
                                        MonomialType const rhs)
    {
        degree_type const degl = degree(lhs);
        degree_type const degr = degree(rhs);

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

    bool is_degree_order() const final { return true; }

    bool is_reverse_order() const final { return true; }

    bool is_weighted_order() const final { return true; }

    params::order type() const final { return params::order::grevlexw; }

    /* vector holding the (positive) weights */
    static std::vector<monomial_base::exponent_type> w_;
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

template <class, class = void>
struct is_weighted_order : std::false_type
{};

/* specialization recognizes monomial order that have nested ::weighted_order */
template <class T>
struct is_weighted_order<T, std::void_t<typename T::weighted_order>>
        : std::true_type
{};

template <class T>
inline constexpr bool is_weighted_order_v = is_weighted_order<T>::value;

}  // namespace gamba
