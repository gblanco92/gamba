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

#include <format>
#include <type_traits>

#include "config.hpp"
#include "container.hpp"
#include "divmask.hpp"
#include "monomial.hpp"
#include "params.hpp"
#include "stats.hpp"

namespace gamba
{

template <class MonomialOrder>
struct spair_type
{
    using monomial_order = MonomialOrder;
    using monomial_type  = monomial<spair_hashtable, monomial_order>;
    using index_type     = monomial_type::index_type;
    using degree_type    = monomial_type::degree_type;

    /* signed type so negative degrees are sorted first */
    using signed_degree_type = std::make_signed_t<degree_type>;

    monomial_type lcm{};
    index_type idx1{};
    index_type idx2{};
    signed_degree_type deg{};
    divmask_type sdm{};
};

template <class MonomialOrder>
class spair_set
{
public:
    using monomial_order = MonomialOrder;
    using spair_type_t   = spair_type<monomial_order>;

    using monomial_context = spair_hashtable;
    using monomial_type    = monomial<monomial_context, monomial_order>;
    using monomial_init    = typename monomial_type::monomial_init_t;

    using basis_monomial_context = basis_hashtable;
    using basis_monomial_type =
        monomial<basis_monomial_context, monomial_order>;

    using const_iterator = aligned_vector<spair_type_t>::const_iterator;

    spair_set()
    {
        queue.reserve(SPAIRSET_INIT_SIZE);
        m_mon_set.reserve(HASHTABLE_INIT_SIZE);
    }

    /* make this class non-copyable */
    spair_set(spair_set const&)            = delete;
    spair_set& operator=(spair_set const&) = delete;

    monomial_type get_lcm(basis_monomial_type const lhs,
                          basis_monomial_type const rhs)
    {
        monomial_init const lcm = monomial_init::lcm(lhs, rhs);

        auto const [it, flag] = m_mon_set.insert(lcm);
        return *it;
    }

    template <class BasisType>
    void rehash_hash_table(BasisType const& basis)
    {
        m_mon_set.clear();

        for (spair_type_t& sp : queue)
        {
            auto const idx1 = static_cast<size_t>(sp.idx1);
            auto const idx2 = static_cast<size_t>(sp.idx2);

            basis_monomial_type const lm1 = basis.lead_mon(idx1);
            basis_monomial_type const lm2 = basis.lead_mon(idx2);

            sp.lcm = get_lcm(lm1, lm2);
        }
    }

    template <class BasisType>
    void update_divmasks(BasisType const& basis)
    {
        /* since the current spairs depend on the generators prior the current
         * update, only update lcm divmasks if divmaks have been updated */
        if (basis.divmasks_version() == m_divmasks_version)
            return;

        /* update & pack divisibility masks for lcm monomials in spairs */
        for (spair_type_t& sp : queue)
        {
            auto const lhs_idx = static_cast<size_t>(sp.idx1);
            auto const rhs_idx = static_cast<size_t>(sp.idx2);

            divmask_type const lhs_sdm = basis.lead_sdm(lhs_idx);
            divmask_type const rhs_sdm = basis.lead_sdm(rhs_idx);

            /* the divmask of a lcm is the bitwise OR of the divmasks */
            sp.sdm.mask = lhs_sdm.mask | rhs_sdm.mask;
        }

        m_divmasks_version = basis.divmasks_version();
    }

    double memory_usage() const
    {
        double mem_size = 0.0;

        mem_size += memory_size(queue);
        mem_size += m_mon_set.memory_usage();

        return mem_size;
    }

    /* align spairs vector to cache lines size to avoid false sharing during
     * multithreaded update process */
    aligned_vector<spair_type_t> queue{};

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
    using spair_type_t      = spair_type<monomial_order>;
    using index_type        = spair_type_t::index_type;
    using signed_index_type = std::make_signed_t<index_type>;

    /* comparators return a signed type so negation for reverse sorting works */
    constexpr int32_t operator()(spair_type_t const& lhs,
                                 spair_type_t const& rhs) const
    {
        if (lhs.lcm != rhs.lcm)
            return monomial_order{}(lhs.lcm, rhs.lcm);

        auto const lhs_idx1 = static_cast<signed_index_type>(lhs.idx1);
        auto const rhs_idx1 = static_cast<signed_index_type>(rhs.idx1);

        /* old indices in new spairs are always in the first component */
        return lhs.deg != rhs.deg ? lhs.deg - rhs.deg : lhs_idx1 - rhs_idx1;
    }
};

template <class MonomialOrder>
struct spair_select_order
{
    using monomial_order = MonomialOrder;
    using spair_type_t   = spair_type<monomial_order>;

    /* comparators return a signed type so negation for reverse sorting works */
    int32_t operator()(spair_type_t const& lhs, spair_type_t const& rhs) const
    {
        if (lhs.deg != rhs.deg)
            return lhs.deg - rhs.deg;

        /* in case of ties use the monomial order */
        return monomial_order{}(lhs.lcm, rhs.lcm);
    }
};

template <class MonomialOrder>
[[nodiscard]] auto select_spairs(spair_set<MonomialOrder>& spairs,
                                 gamba_params const& params)
{
    using monomial_order     = MonomialOrder;
    using spair_type         = spair_type<monomial_order>;
    using degree_type        = spair_type::degree_type;
    using signed_degree_type = spair_type::signed_degree_type;

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    if (params.all_spairs)
    {
        /* if using all spairs just sort by the given monomial order of lcms */
        std::sort(std::begin(spairs.queue), std::end(spairs.queue),
                  [](spair_type const& lhs, spair_type const& rhs) {
                      return monomial_order{}(lhs.lcm, rhs.lcm) < 0;
                  });

        /* timings */
        auto const end_cputime  = std::clock();
        auto const end_walltime = std::chrono::system_clock::now();

        stats.select_walltime +=
            std::chrono::duration<double>(end_walltime - start_walltime)
                .count();
        stats.select_cputime +=
            static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

        stats.spairs_reduced += static_cast<ssize_t>(spairs.queue.size());

        return std::make_tuple(std::ranges::subrange(std::cbegin(spairs.queue),
                                                     std::cend(spairs.queue)),
                               std::numeric_limits<degree_type>::max());
    }

    /* sort first by degree and within same degrees use the monomial order */
    std::sort(std::begin(spairs.queue), std::end(spairs.queue),
              [](spair_type const& lhs, spair_type const& rhs) {
                  return spair_select_order<monomial_order>{}(lhs, rhs) < 0;
              });

    signed_degree_type const min_deg{spairs.queue[0].deg};

    auto const mindeg_end = std::ranges::upper_bound(
        spairs.queue, min_deg, {}, [](spair_type const& sp) { return sp.deg; });

    ssize_t const num_selec_pairs =
        std::min(params.max_spairs, mindeg_end - std::cbegin(spairs.queue));

    auto const spair_end = std::cbegin(spairs.queue) + num_selec_pairs;

    if (params.rounds_info)
    {
        std::cout << std::format("{:3}{:>7} / {:<5}", min_deg, num_selec_pairs,
                                 spairs.queue.size())
                  << std::flush;
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.select_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.select_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    stats.spairs_reduced += std::distance(std::cbegin(spairs.queue), spair_end);

    return std::make_tuple(
        std::ranges::subrange(std::cbegin(spairs.queue), spair_end),
        static_cast<degree_type>(min_deg));
}

}  // namespace gamba
