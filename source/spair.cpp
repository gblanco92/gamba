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

#include "order.hpp"

namespace gamba
{

#ifndef GAMBA_UNITY_BUILD
template class std::pair<size_t, size_t> select_spairs<order_grevlex>(
    spair_set& spairs);
// template class std::pair<size_t, size_t> select_spairs<order_lexic>(
//     spair_set& spairs);
template class std::pair<size_t, size_t> select_spairs<order_blockelim>(
    spair_set& spairs);
#endif

spair_set::spair_set()
{
    queue.reserve(SPAIRSET_INIT_SIZE);
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);
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
