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

#include <flint/flint.h>

#include "base_basis.hpp"
#include "container.hpp"
#include "divmask.hpp"
#include "learn.hpp"
#include "monomial.hpp"
#include "order.hpp"

namespace gamba
{

struct spair_type
{
    using monomial_type = monomial<spair_hashtable>;
    using index_type    = monomial_type::index_type;
    using degree_type   = monomial_type::degree_type;

    /* signed type so negative degrees are sorted first */
    using signed_degree_type = std::make_signed_t<degree_type>;

    monomial_type lcm{};
    index_type idx1{};
    index_type idx2{};
    signed_degree_type deg{};
    divmask_type sdm{};
};

class spair_set
{
public:
    using monomial_context = spair_hashtable;
    using monomial_type    = monomial<monomial_context>;
    using monomial_init    = monomial_type::monomial_init_t;

    using basis_monomial_context = basis_hashtable;
    using basis_monomial_type    = monomial<basis_monomial_context>;

    using const_iterator = aligned_vector<spair_type>::const_iterator;

    spair_set();

    /* make this class non-copyable */
    spair_set(spair_set const&)            = delete;
    spair_set& operator=(spair_set const&) = delete;

    FORCE_INLINE monomial_type get_lcm(basis_monomial_type const lhs,
                                       basis_monomial_type const rhs)
    {
        monomial_init const lcm = monomial_init::lcm(lhs, rhs);

        auto const [it, flag] = m_mon_set.insert(lcm);
        return *it;
    }

    void rehash_hash_table(base_polynomial_basis const& basis);

    void update_divmasks(base_polynomial_basis const& basis);

    double memory_usage() const;

    /* align spairs vector to cache lines size to avoid false sharing during
     * multithreaded update process */
    aligned_vector<spair_type> queue{};

private:
    /* unordered (flat) set storing all monomials appearing as lcm in spairs */
    monomial_set<monomial_type> m_mon_set{};

    /* version of divmap used in the computation of current lcm divmasks */
    size_t m_divmasks_version{0UL};
};

template <class MonomialOrder>
struct spair_update_order
{
    using monomial_order    = MonomialOrder;
    using index_type        = spair_type::index_type;
    using signed_index_type = std::make_signed_t<index_type>;

    /* comparators return a signed type so negation for reverse sorting works */
    constexpr int32_t operator()(spair_type const& lhs,
                                 spair_type const& rhs) const
    {
        if (lhs.lcm != rhs.lcm)
            return monomial_order{}(lhs.lcm, rhs.lcm);

        auto const lhs_idx1 = static_cast<signed_index_type>(lhs.idx1);
        auto const rhs_idx1 = static_cast<signed_index_type>(rhs.idx1);

        /* old indices in new spairs are always in the first component */
        return lhs.deg != rhs.deg ? lhs.deg - rhs.deg : lhs_idx1 - rhs_idx1;
    }
};

struct spair_select_order
{
    /* comparators return a signed type so negation for reverse sorting works */
    int32_t operator()(monomial_order const& mon_order,
                       spair_type const& lhs,
                       spair_type const& rhs) const
    {
        if (lhs.deg != rhs.deg)
            return lhs.deg - rhs.deg;

        /* in case of ties use the monomial order */
        return mon_order.cmp(lhs.lcm, rhs.lcm);
    }
};

[[nodiscard]] std::pair<size_t, size_t> select_spairs(
    spair_set& spairs,
    monomial_order const& mon_order,
    learn_f4_data const* learn_data,
    size_t const round);

}  // namespace gamba
