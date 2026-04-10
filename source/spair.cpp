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

#include "spair.hpp"

#include "learn.hpp"
#include "order.hpp"

namespace gamba
{

spair_set::spair_set()
{
    queue.reserve(SPAIRSET_INIT_SIZE);
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);
}

[[nodiscard]] std::pair<size_t, size_t> select_spairs(
    spair_set& spairs,
    monomial_order const& mon_order,
    learn_f4_data const* learn_data,
    size_t const round)
{
    using degree_type        = spair_type::degree_type;
    using signed_degree_type = spair_type::signed_degree_type;

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    if (params::all_spairs)
    {
        /* if using all spairs just sort by the given monomial order of lcms */
        std::sort(std::begin(spairs.queue), std::end(spairs.queue),
                  [&mon_order](spair_type const& lhs, spair_type const& rhs) {
                      return mon_order.cmp(lhs.lcm, rhs.lcm) < 0;
                  });

        /* timings */
        auto const end_cputime  = std::clock();
        auto const end_walltime = std::chrono::system_clock::now();

        stats::select_walltime +=
            std::chrono::duration<double>(end_walltime - start_walltime)
                .count();
        stats::select_cputime +=
            static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

        stats::spairs_reduced += static_cast<ssize_t>(spairs.queue.size());

        return std::make_pair(spairs.queue.size(),
                              std::numeric_limits<size_t>::max());
    }

    /* sort first by degree and within same degrees use the monomial order */
    std::sort(std::begin(spairs.queue), std::end(spairs.queue),
              [&mon_order](spair_type const& lhs, spair_type const& rhs) {
                  return spair_select_order{}(mon_order, lhs, rhs) < 0;
              });

    signed_degree_type const min_deg{spairs.queue[0].deg};

    auto const mindeg_end = std::ranges::upper_bound(
        spairs.queue, min_deg, {}, [](spair_type const& sp) { return sp.deg; });

    ssize_t /* const */ num_selec_pairs =
        std::min(params::max_spairs, mindeg_end - std::cbegin(spairs.queue));

    log::print(log::INFO2, "│ {:3}{:>8} / {:<5}", min_deg, num_selec_pairs,
               spairs.queue.size());
    ::fflush(stdout);

    /* ignore spairs for current F4 round if all reduce to zero */
    if (learn_data and not learn_data->save_data
        and learn_data->lm_round[round].empty())
    {
        /* remove inserted spairs from spair queue */
        spairs.queue.erase(std::begin(spairs.queue),
                           std::begin(spairs.queue) + num_selec_pairs);

        num_selec_pairs = 0;
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::select_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::select_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    stats::spairs_reduced += num_selec_pairs;

    return std::make_pair(num_selec_pairs, static_cast<degree_type>(min_deg));
}

void spair_set::rehash_hash_table(base_polynomial_basis const& basis)
{
    m_mon_set.clear();

    for (spair_type& sp : queue)
    {
        auto const idx1 = static_cast<size_t>(sp.idx1);
        auto const idx2 = static_cast<size_t>(sp.idx2);

        basis_monomial_type const lm1 = basis.lead_mon(idx1);
        basis_monomial_type const lm2 = basis.lead_mon(idx2);

        sp.lcm = get_lcm(lm1, lm2);
    }
}

void spair_set::update_divmasks(base_polynomial_basis const& basis)
{
    /* since the current spairs depend on the generators prior the current
     * update, only update lcm divmasks if divmaks have been updated */
    if (basis.divmasks_version() == m_divmasks_version)
        return;

    /* update & pack divisibility masks for lcm monomials in spairs */
    for (spair_type& sp : queue)
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

double spair_set::memory_usage() const
{
    double mem_size = 0.0;

    mem_size += memory_size(queue);
    mem_size += m_mon_set.memory_usage();

    return mem_size;
}

}  // namespace gamba
