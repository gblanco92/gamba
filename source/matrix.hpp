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

#include <algorithm>

#include "base_basis.hpp"
#include "divmask.hpp"
#include "kernel/avx2/find.hpp"
#include "monomial.hpp"
#include "spair.hpp"
#include "utils.hpp"

namespace gamba
{

class matrix_f4
{
public:
    using monomial_context = matrix_hashtable;
    using monomial_type    = monomial<monomial_context>;
    using monomial_init    = monomial_type::monomial_init_t;

    static_assert(std::is_standard_layout_v<monomial_type>);

    using exponent_type = monomial_type::exponent_type;
    using index_type    = monomial_type::index_type;
    using count_type    = monomial_type::count_type;

    using basis_type               = base_polynomial_basis;
    using basis_monomial_context   = basis_type::monomial_context;
    using basis_monomial_type      = basis_type::monomial_type;
    using basis_monomial_init      = basis_monomial_type::monomial_init_t;
    using basis_monomial_vect_type = basis_type::monomial_vect_type;

    using divmap_type         = divmask_map<basis_monomial_type>;
    using spair_monomial_type = spair_set::monomial_type;

    using index_vect_type    = std::span<index_type>;
    using monomial_vect_type = std::vector<monomial_type>;
    using monomial_ptr_type  = monomial_type const*;
    using void_ptr_type      = base_polynomial_basis::void_ptr_type;

    using degree_type = basis_type::degree_type;
    using length_type = basis_type::length_type;

    explicit matrix_f4(basis_type& basis) : m_divmap{basis.divmap()} {}

    /* make this class non-copyable */
    matrix_f4(matrix_f4 const&)            = delete;
    matrix_f4& operator=(matrix_f4 const&) = delete;

    /* insert the first num_spairs in of spairs' queue */
    void insert_spairs(spair_set& spairs,
                       size_t const num_spairs,
                       basis_type& basis);

    /* insert generators from basis as top matrix rows for interreduction */
    void insert_generators_reduce(basis_type const& basis);

    void symbolic_preprocessing(basis_type const& basis);

    template <class MonomialOrder>
    void convert_monomials_to_columns(MonomialOrder /*unused*/);

    double memory_usage() const;

    void clear();

    size_t num_cols() const { return m_mon_set.size(); }

    size_t num_rows() const { return m_top_rows.size() + m_bottom_rows.size(); }

    size_t num_top_rows() const { return m_top_rows.size(); }

    size_t num_bottom_rows() const { return m_bottom_rows.size(); }

    monomial_ptr_type top_monomials(size_t const i) const
    {
        return m_top_rows[i].data();
    }

    monomial_ptr_type bottom_monomials(size_t const i) const
    {
        return m_bottom_rows[i].data();
    }

    void_ptr_type top_coeffs(size_t const i) const { return m_top_coefs[i]; }

    void_ptr_type bottom_coeffs(size_t const i) const
    {
        return m_bottom_coefs[i];
    }

    length_type top_size(size_t const i) const
    {
        return static_cast<length_type>(m_top_rows[i].size());
    }

    length_type bottom_size(size_t const i) const
    {
        return static_cast<length_type>(m_bottom_rows[i].size());
    }

    monomial_vect_type& bottom_rows(size_t const i) { return m_bottom_rows[i]; }

private:
    /* the following functions must be templates because polynomial products are
     * performed with monomials from both the spair & the basis hash tables */
    template <class OtherMonomialContext>
    FORCE_INLINE monomial_vect_type create_multiplied_poly_matrix_row(
        basis_monomial_vect_type const& gen_mons,
        monomial<OtherMonomialContext> const& mon);

    template <class OtherMonomialContext>
    FORCE_INLINE monomial_vect_type create_multiplied_poly_matrix_row_prefetch(
        basis_monomial_vect_type const& gen_mons,
        monomial<OtherMonomialContext> const& mon);

    FORCE_INLINE monomial_vect_type
    create_poly_matrix_row(basis_monomial_vect_type const& gen_mons);

    FORCE_INLINE size_t
    find_multiplied_reducer(monomial_type const mon,
                            std::vector<basis_monomial_type> const& gens_lm,
                            aligned_vector<divmask_type> const& gens_sdm) const;

    /* class to generate the divisibility masks (borrowed from basis class) */
    divmap_type const& m_divmap;

    /* sparse top rows of the matrix, the reducers, stored as hashed monomial */
    std::vector<monomial_vect_type> m_top_rows{};

    /* sparse bottom rows of matrix, the reductees, stored as hashed monomial */
    std::vector<monomial_vect_type> m_bottom_rows{};

    /* coefficients of the generators on top rows; memory owned by basis class*/
    std::vector<void_ptr_type> m_top_coefs{};

    /* coefficients of the generators on the bottom rows; mem. owned by basis */
    std::vector<void_ptr_type> m_bottom_coefs{};

    /* unordered (flat) set storing all monomials appearing in the matrix */
    monomial_set<monomial_type> m_mon_set{};

public:
    /* the i-th element of the vector corresponds to the i-th column monomial */
    std::vector<monomial_type> m_col_to_mon{};

    /* new non-zero rows found after reducing bottom rows */
    std::vector<index_vect_type> m_new_rows{};

    /* new non-zero coefs found after reducing bottom rows */
    std::vector<void_ptr_type> m_new_coefs{};
};

/* creates a matrix row by multipliying the monomials in 'gen_mons' by the
 * quotient between 'mon' and the leading monomial of 'gen_mons' */
template <class OtherMonomialContext>
matrix_f4::monomial_vect_type FORCE_INLINE
matrix_f4::create_multiplied_poly_matrix_row(
    basis_monomial_vect_type const& gen_mons,
    monomial<OtherMonomialContext> const& mon)
{
    basis_monomial_type const lead_mon = gen_mons[0];

    /* creates temporary (never inserted in a hash table) monomial quotient
     * in the *basis* hash table buffer since it won't be modified */
    auto const quot = basis_monomial_init::quotient(mon, lead_mon);

    monomial_vect_type row(gen_mons.size());

    /* transform the monomial in polynomial by multiplying them by the quotient
     * monomial and insert the new product monomial in the matrix row */
    std::transform(std::cbegin(gen_mons), std::cend(gen_mons), std::begin(row),
                   [quot, this](basis_monomial_type const gen_mon) {
                       monomial_init const prod_mon =
                           monomial_init::product(quot, gen_mon);

                       /* insert in *matrix* hash table, save hashed index */
                       auto const [it, _] = m_mon_set.insert(prod_mon);

                       return *it;
                   });

    /* mark the lm of the new reducer as a pivot */
    row[0].data().idx = 1;

    return row;
}

template <class OtherMonomialContext>
matrix_f4::monomial_vect_type FORCE_INLINE
matrix_f4::create_multiplied_poly_matrix_row_prefetch(
    basis_monomial_vect_type const& gen_mons,
    monomial<OtherMonomialContext> const& mon)
{
    /* prefetch distance */
    constexpr size_t const PREFETCH_DIST = 2;

    basis_monomial_type const lead_mon = gen_mons[0];

    /* creates temporary (never inserted in a hash table) monomial quotient
     * in the *basis* hash table buffer since it won't be modified */
    auto const quot = basis_monomial_init::quotient(mon, lead_mon);

    monomial_vect_type row(gen_mons.size());

    /* lambda function to multiply and insert basis monomial */
    auto const mult_insert = [quot, this](basis_monomial_type const bmon) {
        monomial_init const prod_mon = monomial_init::product(quot, bmon);

        /* insert in *matrix* hash table, save hashed index */
        auto const [it, _] = m_mon_set.insert(prod_mon);

        return *it;
    };

    auto it        = std::cbegin(gen_mons);
    auto out       = std::begin(row);
    auto const end = std::cend(gen_mons) - PREFETCH_DIST;

    /* edge case for prefetching distance == 2 */
    if (UNLIKELY(gen_mons.size() == 1))
        goto single_monomial;  // NOLINT

    for (; it != end; ++it)
    {
        /* compute product hash at distance PREFETCH_DIST */
        auto const next_it     = it + PREFETCH_DIST;
        size_t const next_hash = next_it->hash() + quot.hash();

        /* use future hash to prefetch first hash table look up */
        m_mon_set.prefetch_insert(next_hash);

        *out++ = mult_insert(*it);
    }

    *out++ = mult_insert(*it++);
single_monomial:
    *out = mult_insert(*it);

    /* mark the lm of the new reducer as a pivot */
    row[0].data().idx = 1;

    return row;
}

FORCE_INLINE matrix_f4::monomial_vect_type matrix_f4::create_poly_matrix_row(
    basis_monomial_vect_type const& gen_mons)
{
    monomial_vect_type row(gen_mons.size());

    /* transform the monomial in polynomial by multiplying them by the quotient
     * monomial and insert the new product monomial in the matrix row */
    std::transform(std::cbegin(gen_mons), std::cend(gen_mons), std::begin(row),
                   [this](basis_monomial_type const gen_mon) {
                       monomial_init const mon = monomial_init::copy(gen_mon);

                       /* insert in *matrix* hash table, save hashed index */
                       auto const [it, _] = m_mon_set.insert(mon);

                       return *it;
                   });

    return row;
}

FORCE_INLINE size_t matrix_f4::find_multiplied_reducer(
    monomial_type const mon,
    std::vector<basis_monomial_type> const& gens_lm,
    aligned_vector<divmask_type> const& gens_sdm) const
{
    /* each matrix monomial' divmask is only computed once, here */
    divmask_type const mon_sdm = m_divmap.compute_divmask(mon);

#if 0
    auto const search_range = std::views::zip(gens_lm, gens_sdm);

    auto const it = std::ranges::find_if(
        search_range,
        /* is the matrix (mon) divisible by a basis leading monomial (lm)? */
        [mon, mon_sdm](auto const lm_sdm) {
            return is_divisible_mask(std::get<0>(lm_sdm),
            std::get<1>(lm_sdm),
                                     mon, mon_sdm);
        });

    return std::distance(std::cbegin(search_range), it);
#endif

    auto const* const lm_masks_ptr =
        reinterpret_cast<divmask_type::mask_type const*>(gens_sdm.data());
    auto const* const lm_ind_ptr =
        reinterpret_cast<monomial_type::index_type const*>(gens_lm.data());

    return find_multiplied_reducer_kernel(
        mon_sdm.mask, mon.cbegin(), monomial_type::exp_size, lm_masks_ptr,
        lm_ind_ptr, basis_monomial_type::exps_vect(), gens_sdm.size());
}

template <class MonomialOrder>
void matrix_f4::convert_monomials_to_columns(MonomialOrder /*unused*/)
{
    using monomial_order = MonomialOrder;

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    /* initialize the column -> monomial mapping with all matrix monomials */
    m_col_to_mon.resize(num_cols());

    std::copy(std::cbegin(m_mon_set), std::cend(m_mon_set),
              std::begin(m_col_to_mon));

    /* sort the monomials appearing as column by *reverse* monomial order */
    std::sort(std::begin(m_col_to_mon), std::end(m_col_to_mon),
              [](monomial_type const lhs, monomial_type const rhs) {
                  return monomial_order{}(lhs, rhs) > 0;
              });

    /* store the inverse mapping: matrix monomial -> column index */
    std::for_each(std::cbegin(m_col_to_mon), std::cend(m_col_to_mon),
                  [idx = 0U](monomial_type const mon) mutable {
                      mon.data().idx = idx++;
                  });

    /* there is no need to sort the individual rows since the rows/generators
     * are sorted by the given monomial order and the stable partition preserves
     * this order */

    /* sort top rows by pivot order */
    std::ranges::sort(
        std::views::zip(m_top_rows, m_top_coefs),
        [](index_type const lhs, index_type const rhs) { return lhs < rhs; },
        [](auto const& row_cfs) { return std::get<0>(row_cfs)[0].data().idx; });

    /* sort bottom rows by (reverse) pivot order */
    std::ranges::sort(
        std::views::zip(m_bottom_rows, m_bottom_coefs),
        [](index_type const lhs, index_type const rhs) { return lhs > rhs; },
        [](auto const& row_cfs) { return std::get<0>(row_cfs)[0].data().idx; });

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::convert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::convert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    size_t const nnz_top = std::accumulate(
        std::cbegin(m_top_rows), std::cend(m_top_rows), 0ULL,
        [](size_t acc, auto const& row) { return acc + row.size(); });

    size_t const nnz_bot = std::accumulate(
        std::cbegin(m_bottom_rows), std::cend(m_bottom_rows), 0ULL,
        [](size_t acc, auto const& row) { return acc + row.size(); });

    double const density =
        static_cast<double>(nnz_top + nnz_bot)
        / static_cast<double>(m_top_rows.size() + m_bottom_rows.size())
        / static_cast<double>(m_mon_set.size()) * 100;

    log::print(log::INFO2, "{:>10} x {:<9} {:>6.2f}%",
               m_top_rows.size() + m_bottom_rows.size(), m_mon_set.size(),
               density);
    ::fflush(stdout);
}

}  // namespace gamba
