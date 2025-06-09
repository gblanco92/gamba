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

#pragma once

#include "flat_set.hpp"

namespace gamba
{

template <class MonomialType>
class divmask_map
{
public:
    using monomial_type = MonomialType;

    using exponent_type     = monomial_type::exponent_type;
    using exponent_ptr_type = exponent_type const*;

    using mask_type = divmask_type::mask_type;

    explicit divmask_map(size_t const num_vars);

    void update_map(monomial_flat_set<monomial_type> const& mon_set);

    template <class AnyMonomialType>
    divmask_type compute_divmask(AnyMonomialType const& mon) const;

    size_t version() const { return m_version_num; }

private:
    static constexpr size_t const m_max_bits = CHAR_BIT * sizeof(divmask_type);

    /* threshold monomial for each bit in m_bits_per_var */
    std::array<size_t, m_max_bits> m_divmap{};

    /* number of variables select for the divmasks */
    size_t m_nvars_divmask{};

    /* number of monomiasl processed in previous update */
    size_t m_prev_size{0};

    /* unique ID for each update of the divmask_map */
    size_t m_version_num{1'000UL};

    /* minimun exponent in variables selects for divmasks */
    std::array<exponent_type, m_max_bits> m_min_deg{};

    /* maximum exponent in variables selects for divmasks */
    std::array<exponent_type, m_max_bits> m_max_deg{};

    /* variables selected for divisibility masks (ignore degree(s)), max = 32 */
    std::array<size_t, m_max_bits> m_vars_divmask{};

    /* number of bits for each selected variable in vars_divmask */
    std::array<size_t, m_max_bits> m_bits_per_var{};
};

template <class MonomialType>
divmask_map<MonomialType>::divmask_map(size_t const num_vars) :
        /* if num_vars > total_bits use the first 'max_bits' variables */
        m_nvars_divmask{std::min(num_vars, m_max_bits)}
{
    std::ranges::fill(m_min_deg, std::numeric_limits<exponent_type>::max());
    std::ranges::fill(m_max_deg, std::numeric_limits<exponent_type>::min());

    /* ignore degree in the first position of exponent vectors */
    size_t idx = 1;

    for (size_t i = 0; i < m_nvars_divmask; ++i)
    {
        /* if non-block order ::block_size is 0 and idx starts at 1 */
        if (idx == monomial_type::block_size)
            ++idx;

        /* if num_vars < total_bits distribute bits across variables with the
         * following simple rule */
        m_bits_per_var[i] = m_max_bits / num_vars + (i < m_max_bits % num_vars);

        /* save variable index */
        m_vars_divmask[i] = idx++;
    }
}

template <class MonomialType>
void divmask_map<MonomialType>::update_map(
    monomial_flat_set<MonomialType> const& mon_set)
{
    /* save a copy of previous min/max exponents to compare with updated */
    std::array<exponent_type, m_max_bits> const prev_min_deg{m_min_deg};
    std::array<exponent_type, m_max_bits> const prev_max_deg{m_max_deg};

    auto const start = mon_set.cbegin() + m_prev_size;
    auto const end   = mon_set.cend();

    /* update min/max of each exponent across all monomials in the hash table */
    for (auto mon_it = start; mon_it != end; ++mon_it)
    {
        exponent_ptr_type const exp = mon_it->cbegin();

        for (size_t i = 0; i < m_nvars_divmask; ++i)
        {
            m_min_deg[i] = std::min(m_min_deg[i], exp[m_vars_divmask[i]]);

            m_max_deg[i] = std::max(m_max_deg[i], exp[m_vars_divmask[i]]);
        }
    }

    /* save number of processed monomials so far for next call to update_map */
    m_prev_size = mon_set.size();

    /* if no changes in min/max exponents do nothing and return */
    if (std::ranges::equal(m_min_deg, prev_min_deg)
        and std::ranges::equal(m_max_deg, prev_max_deg))
    {
        return;
    }

    /* update cut-points; divide the range [min, max] into bits_per_var + 1
     * equal pieces and use the #bits_per_var midpoints of those ranges as the
     * cut points for the bits */
    for (size_t i = 0, j = 0; i < m_nvars_divmask; ++i)
    {
        exponent_type inc =
            static_cast<exponent_type>(m_max_deg[i] - m_min_deg[i] + 1)
            / static_cast<exponent_type>(m_bits_per_var[i] + 1);

        if (inc == 0)
            inc = 1;

        for (size_t k = 1; k <= m_bits_per_var[i]; ++k)
        {
            m_divmap[j++] = static_cast<exponent_type>(m_min_deg[i] + inc * k);
        }
    }

    /* update divmask_map version number */
    m_version_num++;
}

template <class MonomialType>
template <class OtherMonomialType>
[[nodiscard]] divmask_type divmask_map<MonomialType>::compute_divmask(
    OtherMonomialType const& mon) const
{
    exponent_ptr_type const exps = mon.cbegin();

    divmask_type sdm{0};

    for (size_t i = 0, bit = 0; i < m_nvars_divmask; ++i)
    {
        /* get exponent of the i-th select variable for divmask computation */
        exponent_type const ei = exps[m_vars_divmask[i]];

        /* if max = min + 1, then inc = 1 and divmap[] = max, therefore use >=
         * instead of > to avoid a constant useless bit in the all divmasks */
        for (size_t k = 0; k < m_bits_per_var[i]; ++k)
        {
            // sdm.mask |= (ei >= m_divmap[bit]) << bit++;
            sdm.mask = (sdm.mask << 1U) | (ei >= m_divmap[bit++]);
        }
    }

    return sdm;
}

}  // namespace gamba
