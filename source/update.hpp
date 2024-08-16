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

#include <limits>

#include "divmask.hpp"
#include "spair.hpp"
#include "stats.hpp"

namespace gamba
{

template <class MonomialOrder, class BasisType>
void update_spairs(spair_set<MonomialOrder>& spairs,
                   BasisType& basis,
                   size_t const h) /* index of an element in basis */
{
    using monomial_order = MonomialOrder;
    using basis_type     = BasisType;
    using spair_set_type = ::gamba::spair_set<monomial_order>;

    using monomial_context = spair_set_type::monomial_context;
    using monomial_type    = spair_set_type::monomial_type;

    using basis_monomial_context = basis_type::monomial_context;
    using basis_monomial_type    = basis_type::monomial_type;

    using spair_type         = spair_type<monomial_order>;
    using degree_type        = spair_type::degree_type;
    using index_type         = basis_type::index_type;
    using signed_degree_type = spair_type::signed_degree_type;

    static_assert(std::is_signed_v<signed_degree_type>);

    static constexpr auto is_divisible =
        monomial_type::template check_monomial_division_divmask<
            monomial_context>;

    static constexpr auto is_divisible_other =
        monomial_type::template check_monomial_division_divmask<
            basis_monomial_context>;

    /* we assume that the basis elements with indices [0, h) have
     * already been used to create spairs previously & the generator with index
     * h is not deleted */
    assert(h < basis.num_gens());
    assert(not basis.is_deleted(h));

    basis_monomial_type const lm_h = basis.lead_mon(h);
    divmask_type const sdm_h       = basis.lead_sdm(h);

    size_t const prev_size = spairs.queue.size();

    spairs.queue.reserve(prev_size + h);

    // TODO parallel: to avoid false sharing: chunk size == cache line size
    // TODO parallel: each iteration of this loop is independent; chunk size
    /* create new spairs involving the new element */
    for (size_t i = 0; i < h; ++i)
    {
        /* no need to generate spairs for deleted generators in basis */
        if (basis.is_deleted(i))
            continue;

        basis_monomial_type const lm_i = basis.lead_mon(i);
        divmask_type const sdm_i       = basis.lead_sdm(i);

        /* generate lcm of leading terms inside basis hashtable */
        monomial_type const lcm = spairs.get_lcm(lm_i, lm_h);
        /* the divmask of a lcm is the bitwise OR of the divmasks */
        divmask_type const sdm = divmask_type{sdm_i.mask | sdm_h.mask};
        /* the degree of the least common multiple monomial */
        degree_type const deg_lcm = lcm.degree();

        /* create the new spair */
        spairs.queue.push_back({.lcm  = lcm,
                                .idx1 = static_cast<index_type>(i),
                                .idx2 = static_cast<index_type>(h),
                                .deg = static_cast<signed_degree_type>(deg_lcm),
                                .sdm = sdm});
        spair_type& sp_i = spairs.queue.back();

        /* increase spair count for each generator in spair */
        basis.spair_count(i)++;
        basis.spair_count(h)++;

        /* even if the element is redundant we need to create a spair since old
         * spairs may use this redundant element and we needed to have its lcm
         * for Gebauer-Moeller old spairs handling */
        if (basis.is_redundant(i))
        {
            sp_i.deg = -1;

            continue;
        }

        /* set degree to -2 so it is placed at the beginning of all spairs with
         * the same lcm when sorting the new spairs */
        if (basis_monomial_type::coprime_monomials(lm_i, lm_h))
        {
            sp_i.deg = -2;

            continue;
        }

        /* compute spair degree, faster if the order is a degree order */
        if constexpr (is_degree_order_v<monomial_order>)
        {
            sp_i.deg = static_cast<signed_degree_type>(deg_lcm);
        }
        else
        {
            auto const deg_i = static_cast<signed_degree_type>(
                (deg_lcm - lm_i.degree()) + basis.degree(i));
            auto const deg_h = static_cast<signed_degree_type>(
                (deg_lcm - lm_h.degree()) + basis.degree(h));

            sp_i.deg = deg_i > deg_h ? deg_i : deg_h;
        }

        assert(sp_i.deg >= 0);
    }

    auto const new_sp =
        std::begin(spairs.queue) + static_cast<ssize_t>(prev_size);

    // TODO parallel: each iteration of this loop is independent; chunk size
    /* Gebauer-Moeller: check old spairs first */
    for (size_t i = 0; i < prev_size; ++i)
    {
        /* there can't be old spairs not yet deleted */
        assert(spairs.queue[i].deg != -1);

        /* here we are using non-deleted polynomial indices to address the
         * correct non-deleted spair in new_sp range */
        index_type const lhs_del_idx =
            basis.nondel_gen_idx(spairs.queue[i].idx1);
        index_type const rhs_del_idx =
            basis.nondel_gen_idx(spairs.queue[i].idx2);

        /* generators defining an spair cannot be deleted */
        assert(lhs_del_idx != std::numeric_limits<index_type>::max());
        assert(rhs_del_idx != std::numeric_limits<index_type>::max());

        monomial_type const lcm_rl = spairs.queue[i].lcm;
        divmask_type const sdm_rl  = spairs.queue[i].sdm;

        assert(lhs_del_idx < spairs.queue.size() - prev_size);
        assert(rhs_del_idx < spairs.queue.size() - prev_size);

        /* lcm(lhs, h) != lcm(lhs, rhs) and lcm(rhs, h) != lcm(lhs, rhs) */
        /* and lm(h) | lcm(lhs, rhs) */
        if (new_sp[lhs_del_idx].lcm != lcm_rl
            and new_sp[rhs_del_idx].lcm != lcm_rl
            and is_divisible_other(lm_h, sdm_h, lcm_rl, sdm_rl))
        {
            spairs.queue[i].deg = -1;
        }
    }

    // TODO parallel: policy
    /* sort new spairs by increasing order of lcm monomials, if spairs have the
     * same lcm sort them by smallest degree (coprime monomials coming first),
     * finally, in case of a tie prefer earlier polynomials */
    std::sort(new_sp, std::end(spairs.queue),
              [](spair_type const& lhs, spair_type const& rhs) {
                  return spair_update_order<monomial_order>{}(lhs, rhs) < 0;
              });

    // TODO parallel: each iteration of this loop is independent; even if the
    // j-th spair that invalidates the i-th spair is deleted first by another
    // thread, the k-th spair that invalidates the j-th pair must also
    // invalidate the i-th pair (by transitivity of divisibility); chunk size
    /* Gebauer-Moeller: check strict multiples of news pairs */
    for (size_t i = prev_size; i < spairs.queue.size(); ++i)
    {
        /* don't eliminate a coprime spair at this stage, it will be eliminated
         * at the end anyway and if deleted now we can't use it at the next step
         * and we might miss some spair deletions when lcms are equal */
        if (spairs.queue[i].deg < 0)
            continue;

        monomial_type const lcm_i = spairs.queue[i].lcm;
        divmask_type const sdm_i  = spairs.queue[i].sdm;

        /* since a monomial order refines the partial order of divisibility and
         * the spairs are sorted by lcm, it is enough to inspect the smaller new
         * spairs lcm(j, h) | lcm(i, h) (strictly) ==> lcm(j, h) < lcm(i, h) */
        for (size_t j = prev_size; j < i; ++j)
        {
            if (spairs.queue[j].deg == -1)
                continue;

            monomial_type const lcm_j = spairs.queue[j].lcm;
            divmask_type const sdm_j  = spairs.queue[j].sdm;

            /* lcm(i, h) != lcm(j, h) and lcm(j, h) | lcm(i, h) */
            if (lcm_i != lcm_j and is_divisible(lcm_j, sdm_j, lcm_i, sdm_i))
            {
                spairs.queue[i].deg = -1;

                break;
            }
        }
    }

    // TODO parallel: each loop is independent; the first (not yet deleted)
    // spair in a series of spairs with same lcm will never be deleted, and this
    // first spair is the one that decides if the rest of spairs with same lcm
    // are deleted; chunk size
    /* Gebauer-Moeller: check new spairs with equal lcm */
    for (size_t i = prev_size + 1; i < spairs.queue.size(); ++i)
    {
        /* coprime spairs will be eliminated anyway so don't bother */
        if (spairs.queue[i].deg < 0)
            continue;

        monomial_type const lcm_i = spairs.queue[i].lcm;

        /* delete any spair that is not the first in a series of new spairs with
         * the same lcm; such a first element might be a coprime spair and since
         * coprime spairs are always placed first, product criterion applies */
        for (size_t j = i; j > prev_size and spairs.queue[j - 1].lcm == lcm_i;
             --j)
        {
            if (spairs.queue[j - 1].deg == -1)
                continue;

            spairs.queue[i].deg = -1;

            break;
        }
    }

    /* move redundant spairs at the end of spairs queue (do not use
     * std::remove_if if need to access to redundant spairs )*/
    auto const new_end =
        std::partition(std::begin(spairs.queue), std::end(spairs.queue),
                       [](spair_type const& sp) { return sp.deg >= 0; });

    /* decrease spair count for each generator in a redundant spair */
    std::for_each(new_end, std::end(spairs.queue),
                  [&basis](spair_type const& sp) {
                      basis.spair_count(sp.idx1)--;
                      basis.spair_count(sp.idx2)--;
                  });

    /* every time spair_count are decreased new deleted gens. can be created */
    basis.update_deleted_gens();

    stats.gm_criteria += std::distance(new_end, std::end(spairs.queue));

    spairs.queue.erase(new_end, std::cend(spairs.queue));
}

template <class MonomialOrder, class BasisType>
void update_f4(spair_set<MonomialOrder>& spairs,
               BasisType& basis,
               size_t const prev_num_gens)
{
    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    /* update divmasks in spairs */
    spairs.update_divmasks(basis);

    /* update and remove spairs */
    for (size_t i = prev_num_gens; i < basis.num_gens(); ++i)
    {
        update_spairs(spairs, basis, i);
    }

    /* trigger a rehash of spair's hash table */
    if (prev_num_gens != basis.num_gens())
    {
        spairs.rehash_hash_table(basis);
    }

    /* update and make generators redundant */
    basis.update_reduced_gens(prev_num_gens);

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.update_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.update_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

}  // namespace gamba
