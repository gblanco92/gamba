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

#include <vector>

#include "container.hpp"
#include "divmap.hpp"
#include "monomial.hpp"
#include "order.hpp"

namespace gamba
{

class base_polynomial_basis
{
public:
    using monomial_context = basis_hashtable;
    using monomial_type    = monomial<monomial_context>;
    using monomial_init    = monomial_type::monomial_init_t;

    using hash_type      = monomial_type::hash_type;
    using degree_type    = monomial_type::degree_type;
    using index_type     = monomial_type::index_type;
    using index_ptr_type = index_type const*;
    using void_ptr_type  = void const*;

    using length_type = uint32_t;
    using count_type  = uint32_t;

    using divmap_type = divmask_map<monomial_type>;

    using matrix_monomial_context = matrix_hashtable;
    using matrix_monomial_type    = monomial<matrix_monomial_context>;

    using monomial_vect_type = std::span<monomial_type>;
    using index_vect_type    = std::span<index_type>;

protected:
    static constexpr auto is_divisible =
        monomial_type::template check_monomial_division_divmask<
            monomial_context>;

    explicit base_polynomial_basis(size_t const prime) :
            m_field_char{prime},
            m_mon_set{nullptr},
            m_divmap{nullptr}
    {}

    base_polynomial_basis(size_t const n_vars, size_t const f_char);

    void clear();

    void clear_redundant();

    double memory_usage() const;

public:
    virtual ~base_polynomial_basis() = default;

    /* make this class non-copyable */
    base_polynomial_basis(base_polynomial_basis const&)            = delete;
    base_polynomial_basis& operator=(base_polynomial_basis const&) = delete;

    base_polynomial_basis(base_polynomial_basis&&) noexcept = default;
    base_polynomial_basis& operator=(base_polynomial_basis&&) noexcept =
        default;

    /* access member functions */
    virtual uint32_t field_char() const = 0;

    size_t num_gens() const { return m_num_gens; }

    length_type length(size_t const i) const
    {
        return static_cast<length_type>(m_mons[i].size());
    }

    degree_type degree(size_t const i) const { return m_degs[i]; }

    monomial_type lead_mon(size_t const i) const { return m_mons[i][0]; }

    divmask_type lead_sdm(size_t const i) const
    {
        /* lead_sdm must always be computed using latest version of 'divmap' */
        assert(m_divmasks_version == m_divmap->version());

        return m_lead_sdm[i];
    }

    std::vector<index_type> const& reduced_indices() const
    {
        return m_reduced_gens;
    }

    bool is_redundant(size_t const i) const { return m_redundant[i]; }

    bool is_deleted(size_t const i) const
    {
        return m_redundant[i] and m_spair_count[i] == 0;
    }

    size_t num_nondel_gens() const { return m_num_nondel_gens; }

    index_type nondel_gen_idx(size_t const i) const { return m_nondel_gens[i]; }

    bool is_trivial() const { return m_is_trivial; }

    divmap_type const& divmap() const { return *m_divmap; }

    monomial_vect_type const& monomials(size_t const i) const
    {
        return m_mons[i];
    }

    virtual void_ptr_type v_coefficients(size_t const i) const = 0;

    size_t divmasks_version() const { return m_divmasks_version; }

    count_type& spair_count(size_t const i) { return m_spair_count[i]; }

    /* after an update process enforce the following invariant:
     * redundant[reduced_gens[i]] == false */
    template <class MonomialOrder>
    void update_reduced_gens(size_t const prev_num_gens,
                             MonomialOrder /*unused*/);

    void remove_redundant_gens(
        [[maybe_unused]] monomial_order const& mon_order);

    void update_deleted_gens();

    virtual void v_free_generator(size_t const i) = 0;

    void print_info() const;

protected:
    /* number of variables in the polynomial ring */
    size_t m_num_vars{};

    /* field characteristic */
    size_t m_field_char{};

    /* numer of total generators (non-reduced) in the basis */
    size_t m_num_gens{};

    /* unordered (flat) set storing all monomials appearing in the basis; this
     * member variable is *shared* between different bases */
    std::shared_ptr<monomial_set<monomial_type>> m_mon_set{};

    /* hashed monomials for each generator in the current basis */
    std::vector<monomial_vect_type> m_mons{};

    /* total degree of each polynomial in the basis */
    std::vector<degree_type> m_degs{};

    /* helper class to generate the divisibility masks; this member variable is
     * *shared* between different bases */
    std::shared_ptr<divmap_type> m_divmap;

    /* divisibility mask for each generator's leading monomial */
    mutable std::vector<divmask_type> m_lead_sdm{};

    /* version of divmap used in the computation of current lead_sdm */
    mutable size_t m_divmasks_version{0UL};

    /* polynomials of the basis made redundant by G-M update; make it aligned so
     * it can be updated concurrently, avoid std::vector<bool> specialization */
    aligned_vector<uint8_t> m_redundant{};

    /* keeps the index of top reduced polynomials in the basis; this means that
     * for all index i: redundant[reduced_gens[i]] == false for all i */
    std::vector<index_type> m_reduced_gens{};

    /* number of non-deleted generators in the basis; a generators is marked as
     * deleted if it is redundant and does not appear in an spair */
    size_t m_num_nondel_gens{};

    /* number of spairs where the i-th polynomial of the basis appears */
    std::vector<count_type> m_spair_count{};

    /* if a generator is not deleted this maps it to its relative index within
     * the non-deleted generators; otherwise is infinity */
    std::vector<index_type> m_nondel_gens{};

    /* whether the basis is homogeneous or not */
    bool m_is_homogeneous{};

    /* the ideal is equal to the whole ring */
    bool m_is_trivial{};
};

template <class MonomialOrder>
void base_polynomial_basis::update_reduced_gens(size_t const prev_num_gens,
                                                MonomialOrder /*unused*/)
{
    using monomial_order = MonomialOrder;

    size_t num_old_redundant = 0;

    // TODO: parallel, each loop iteration is independent of each other;
    /* check redundancy of old reduced elements in the basis */
    for (index_type const red_idx : m_reduced_gens)
    {
        /* previous reduced elements in basis cannot be already redundant */
        assert(not m_redundant[red_idx]);

        monomial_type const lm_ri = lead_mon(red_idx);
        divmask_type const sdm_ri = lead_sdm(red_idx);

        [[maybe_unused]] degree_type const deg_ri = degree(red_idx);
        [[maybe_unused]] degree_type const deg_lm_ri =
            monomial_order::degree(lm_ri);

        /* if there exist an element from current round dividing the old
         * i-th element, the i-th element can be made redundant */
        for (size_t j = prev_num_gens; j < m_num_gens; ++j)
        {
            monomial_type const lm_j = lead_mon(j);
            divmask_type const sdm_j = lead_sdm(j);

            [[maybe_unused]] degree_type const deg_j = degree(j);
            [[maybe_unused]] degree_type const deg_lm_j =
                monomial_order::degree(lm_j);

            if (is_divisible(lm_j, sdm_j, lm_ri, sdm_ri))
            {
#if NONDEG_ORDERS_HEURISTIC
                if constexpr (not is_degree_order_v<monomial_order>)
                {
                    if (deg_ri - deg_lm_ri < deg_j - deg_lm_j)
                        continue;
                }
#endif
                m_redundant[red_idx] = true;
                num_old_redundant++;

                break;
            }
        }
    }

    // TODO: parallel, each loop iteration is independent of each other
    /* check redundancy of element from the same update batch */
    for (size_t i = prev_num_gens; i < m_num_gens; ++i)
    {
        monomial_type const lm_i = lead_mon(i);
        divmask_type const sdm_i = lead_sdm(i);

        [[maybe_unused]] degree_type const deg_i = degree(i);
        [[maybe_unused]] degree_type const deg_lm_i =
            monomial_order::degree(lm_i);

        /* if there exist an element updated later that divides the current i-th
         * element, the i-th element can be made redundant */
        for (size_t j = i + 1; j < m_num_gens; ++j)
        {
            monomial_type const lm_j = lead_mon(j);
            divmask_type const sdm_j = lead_sdm(j);

            [[maybe_unused]] degree_type const deg_j = degree(j);
            [[maybe_unused]] degree_type const deg_lm_j =
                monomial_order::degree(lm_j);

            if (is_divisible(lm_j, sdm_j, lm_i, sdm_i))
            {
#if NONDEG_ORDERS_HEURISTIC
                if constexpr (not is_degree_order_v<monomial_order>)
                {
                    if (deg_i - deg_lm_i < deg_j - deg_lm_j)
                        continue;
                }
#endif
                m_redundant[i] = true;

                break;
            }
        }
    }

    /* avoid unnecessary copying if there no redundant elements are found */
    if (num_old_redundant > 0)
    {
        /* remove redundant elements from the old reduced part */
        auto const new_end = std::remove_if(
            std::begin(m_reduced_gens), std::end(m_reduced_gens),
            [this](index_type const idx) { return m_redundant[idx]; });

        m_reduced_gens.erase(new_end, std::cend(m_reduced_gens));

        /* reserve enough memory for the new top reduced elements */
        m_reduced_gens.reserve(m_reduced_gens.size() + m_num_gens
                               - prev_num_gens);
    }

    /* store the indices of non-redundant new elements in the basis */
    std::ranges::copy_if(
        std::views::iota(prev_num_gens, m_num_gens),
        std::back_inserter(m_reduced_gens),
        [this](size_t const idx) { return not m_redundant[idx]; });

    stats::redundant_elem =
        static_cast<ssize_t>(m_num_gens - m_reduced_gens.size());
}

}  // namespace gamba
