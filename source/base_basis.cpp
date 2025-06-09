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

#include "base_basis.hpp"

#include "logger.hpp"

namespace gamba
{

base_polynomial_basis::base_polynomial_basis(size_t const n_vars,
                                             size_t const f_char) :
        m_num_vars{n_vars},
        m_field_char{f_char},
        m_mon_set{std::make_shared<monomial_set<monomial_type>>()},
        m_divmap{std::make_shared<divmap_type>(n_vars)}
{}

void base_polynomial_basis::clear()
{
    m_num_gens = 0UL;

    m_mons.clear();
    m_mons.shrink_to_fit();

    m_degs.clear();
    m_degs.shrink_to_fit();

    m_lead_sdm.clear();
    m_lead_sdm.shrink_to_fit();

    m_redundant.clear();
    m_redundant.shrink_to_fit();

    m_reduced_gens.clear();
    m_reduced_gens.shrink_to_fit();

    m_num_nondel_gens = 0UL;

    m_spair_count.clear();
    m_spair_count.shrink_to_fit();

    m_nondel_gens.clear();
    m_nondel_gens.shrink_to_fit();
}

void base_polynomial_basis::remove_redundant_gens()
{
    /* when inserting elements into the basis in decreasing order of leading
     * monomial the Gebauer-Moeller installation ensures that the leading
     * monomials of the final basis are already reduced; if a different order is
     * chosen generators the non-reduced lead mons. must be removed manually */
#if INSERT_ELEMENTS_DECREASING == 0
    constexpr bool const remove_redundant = true;
    /* if the non-redundant heuristic for non-degree orderings is enable AND the
     * order is non-degree we must remove potencial redudant generators */
#elif NONDEG_ORDERS_HEURISTIC == 1
    constexpr bool const remove_redundant =
        not is_degree_order_v<MonomialOrder>;
#else
    constexpr bool const remove_redundant = false;
#endif

    if constexpr (not remove_redundant)
        return;

    size_t num_redundant = 0;

    /* remove redundant elements in final non-reduced Groebner basis */
    for (index_type const i : m_reduced_gens)
    {
        /* reduced elements in basis cannot be already redundant */
        assert(not m_redundant[i]);

        monomial_type const lm_ri = lead_mon(i);
        divmask_type const sdm_ri = lead_sdm(i);

        /* check if there is a reduced element in the basis dividing the i-th
         * reduced element (different from itself) */
        for (index_type const j : m_reduced_gens)
        {
            if (j == i)
                continue;

            if (is_divisible(lead_mon(j), lead_sdm(j), lm_ri, sdm_ri))
            {
                m_redundant[i] = true;
                num_redundant++;

                break;
            }
        }
    }

    if (num_redundant > 0)
    {
        /* do the actual removal of redundant elements indices */
        auto const new_end = std::remove_if(
            std::begin(m_reduced_gens), std::end(m_reduced_gens),
            [this](index_type const idx) { return m_redundant[idx]; });

        m_reduced_gens.erase(new_end, std::cend(m_reduced_gens));
    }

    stats::redundant_elem =
        static_cast<ssize_t>(m_num_gens - m_reduced_gens.size());
}

void base_polynomial_basis::update_deleted_gens()
{
    /* dummy value representing non-deleted generators */
    constexpr index_type const infty = std::numeric_limits<index_type>::max();

    m_num_nondel_gens = 0UL;

    for (size_t i = 0; i < m_num_gens; ++i)
    {
        m_nondel_gens[i] =
            (is_deleted(i) ? infty
                           : static_cast<index_type>(m_num_nondel_gens++));
    }
}

double base_polynomial_basis::memory_usage() const
{
    double mem_size = 0.0;

    mem_size += memory_size(m_mons);
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        mem_size += memory_size(m_mons[i]);
    }

    mem_size += m_mon_set->memory_usage();

    mem_size += memory_size(m_degs);

    mem_size += memory_size(m_redundant);

    mem_size += memory_size(m_reduced_gens);

    mem_size += memory_size(m_spair_count);

    mem_size += memory_size(m_nondel_gens);

    return mem_size;
}

void base_polynomial_basis::print_info() const
{
    static bool printed = false;

    if (printed)
        return;

    size_t const num_mons = std::accumulate(
        std::cbegin(m_mons), std::cend(m_mons), 0ULL,
        [](size_t acc, auto const& v) { return acc + v.size(); });

    auto const max_degree_poly =
        std::max_element(std::cbegin(m_degs), std::cend(m_degs));

    size_t const max_degree = (max_degree_poly != std::cend(m_degs)
                                   ? *max_degree_poly
                                   : std::numeric_limits<degree_type>::max());

    auto const avg_len =
        static_cast<double>(m_num_gens) / static_cast<double>(num_mons);

    log::print(log::INFO1, "\n┌{0:─^{1}}┐\n", " BASIS INFO ", 38);
    log::print(log::INFO1, "│ num. variables: {:>20} │\n",
               fmt::format("{:<10}", m_num_vars));
    log::print(log::INFO1, "│ field charac.: {:>21} │\n",
               fmt::format("{:<10}", m_field_char));
    log::print(log::INFO1, "│ num. generators: {:>19} │\n",
               fmt::format("{:<10}", m_num_gens));
    log::print(log::INFO1, "│ num. monomials: {:>20} │\n",
               fmt::format("{:<10}", num_mons));
    log::print(log::INFO1, "│ avg. length: {:>23} │\n",
               fmt::format("{:<10.2f}", avg_len));
    log::print(log::INFO1, "│ max. degree: {:>23} │\n",
               fmt::format("{:<10}", max_degree));
    log::print(log::INFO1, "│ homogeneous: {:>23} │\n",
               fmt::format("{:<10}", m_is_homogeneous));
    log::print(log::INFO1, "└{0:─^{1}}┘\n", "", 38);

    printed = true;
}

}  // namespace gamba
