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
#include <limits>

#include <immintrin.h>

#include "basis.hpp"
#include "config.hpp"
#include "divmask.hpp"
#include "kernel/find.hpp"
#include "linalg_v3.hpp"
#include "monomial.hpp"
#include "params.hpp"
#include "spair.hpp"

namespace gamba
{

template <class CoefficientType, class MonomialOrder>
class matrix_f4
{
public:
    using coeff_type     = CoefficientType;
    using coeff_type_ptr = coeff_type const*;
    using monomial_order = MonomialOrder;

    using monomial_context = matrix_hashtable;
    using monomial_type    = monomial<monomial_context, monomial_order>;
    using monomial_init    = typename monomial_type::monomial_init_t;

    static_assert(std::is_standard_layout_v<monomial_type>);

    using exponent_type = typename monomial_type::exponent_type;
    using index_type    = typename monomial_type::index_type;
    using count_type    = typename monomial_type::count_type;

    using basis_type             = polynomial_basis<coeff_type, monomial_order>;
    using basis_monomial_context = typename basis_type::monomial_context;
    using basis_monomial_type    = typename basis_type::monomial_type;
    using basis_monomial_init = typename basis_monomial_type::monomial_init_t;
    using basis_monomial_vect_type = typename basis_type::monomial_vect_type;

    using spair_type_t = spair_type<monomial_order>;
    using divmap_type  = divmask_map<monomial_order, basis_monomial_type>;

    using spair_set_type      = spair_set<monomial_order>;
    using spair_monomial_type = spair_set_type::monomial_type;

    using coeff_vect_type    = basis_type::coeff_vect_type;
    using index_vect_type    = std::span<index_type>;
    using monomial_vect_type = std::vector<monomial_type>;
    using monomial_type_ptr  = monomial_type const*;

    using degree_type = basis_type::degree_type;
    using length_type = basis_type::length_type;

private:
    /* we need to check if a matrix monomial is divisible by a basis monomial */
    static constexpr auto divides =
        monomial_type::template check_monomial_division<basis_monomial_context>;

    static constexpr auto divides_mask =
        monomial_type::template check_monomial_division_divmask<
            basis_monomial_context>;

public:
    explicit matrix_f4(basis_type& basis) :
            m_field{basis.field},
            m_divmap{basis.divmap()}
    {}

    /* make this class non-copyable */
    matrix_f4(matrix_f4 const&)            = delete;
    matrix_f4& operator=(matrix_f4 const&) = delete;

    /* insert a range of spairs */
    template <std::ranges::random_access_range RangeType>
    void insert_spairs(spair_set_type& spairs,
                       RangeType spair_range,
                       basis_type& basis);

    void symbolic_preprocessing(basis_type const& basis);

    void convert_monomials_to_columns(gamba_params const& params);

    void convert_monomials_to_columns_reduce(gamba_params const& params);

    void sort_top_rows();

    void sort_bottom_rows();

    size_t num_cols() const { return m_mon_set.size(); }

    size_t num_rows() const { return m_top_rows.size() + m_bottom_rows.size(); }

    size_t num_top_rows() const { return m_top_rows.size(); }

    size_t num_bottom_rows() const { return m_bottom_rows.size(); }

    monomial_type_ptr top_monomials(size_t const i) const
    {
        return m_top_rows[i].data();
    }

    monomial_type_ptr bottom_monomials(size_t const i) const
    {
        return m_bottom_rows[i].data();
    }

    coeff_type_ptr top_coeffs(size_t const i) const { return m_top_coeffs[i]; }

    coeff_type_ptr bottom_coeffs(size_t const i) const
    {
        return m_bottom_coeffs[i];
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

    void insert_generators_reduce(basis_type const& basis);

    void extract_new_rows(linalg_v3<coeff_type, monomial_order>& linalg);

    double memory_usage() const;

    void clear();

private:
    template <class OtherMonomialContext>
    monomial_vect_type create_multiplied_poly_matrix_row(
        basis_monomial_vect_type const& gen_mons,
        monomial<OtherMonomialContext, MonomialOrder> const& mon);

    template <class OtherMonomialContext>
    monomial_vect_type create_multiplied_poly_matrix_row_prefetch(
        basis_monomial_vect_type const& gen_mons,
        monomial<OtherMonomialContext, MonomialOrder> const& mon);

    monomial_vect_type create_poly_matrix_row(
        basis_monomial_vect_type const& gen_mons);

    size_t find_multiplied_reducer(
        monomial_type const mon,
        std::vector<basis_monomial_type> const& gens_lm,
        aligned_vector<divmask_type> const& gens_sdm) const;

    template <class, class>
    friend class polynomial_basis;

    template <class, class>
    friend class linalg_v3;

public:
    field_traits<coeff_type> const m_field;

private:
    /* sparse top rows of the matrix, the reducers, stored as hashed monomial */
    std::vector<monomial_vect_type> m_top_rows{};

    /* sparse bottom rows of matrix, the reductees, stored as hashed monomial */
    std::vector<monomial_vect_type> m_bottom_rows{};

    /* coefficients of the generators on top rows; memory owned by basis class*/
    std::vector<coeff_type_ptr> m_top_coeffs{};

    /* coefficients of the generators on the bottom rows; mem. owned by basis */
    std::vector<coeff_type_ptr> m_bottom_coeffs{};

    /* unordered (flat) set storing all monomials appearing in the matrix */
    monomial_set<monomial_type> m_mon_set{};

    /* the i-th element of the vector corresponds to the i-th column monomial;
     * this is not a row, so no need to align this = monomial_vect_type */
    std::vector<monomial_type> m_col_to_mon{};

    /* class to generate the divisibility masks (borrowed from basis class) */
    divmap_type const& m_divmap;

public:
    /* new non-zero rows found after reducing bottom rows */
    std::vector<index_vect_type> m_new_rows{};

    /* new non-zero coefs found after reducing bottom rows */
    std::vector<coeff_vect_type> m_new_coefs{};
};

template <class CoefficientType, class MonomialOrder>
template <std::ranges::random_access_range RangeType>
void matrix_f4<CoefficientType, MonomialOrder>::insert_spairs(
    spair_set_type& spairs,
    RangeType spair_range,
    polynomial_basis<CoefficientType, MonomialOrder>& basis)
{ /*
   * The matrix insertion of spairs involves the following steps:
   * 1. Allocate memory for hash table and other variables.
   * 2. Sort spairs by increasing monomial order of the lcms.
   * 3. For each range of spairs with given lcms:
   * 3.1. Get all the unique indices of the generators giving these spairs.
   * 3.2. The shortest of such generators goes into the top part of the matrix.
   * 3.3. The other generators go to the lower part of the matrix.
   * 3.4. Save top/bottom generator indices to access coefficients later on.
   * 3.5. Clear memory and remove unused spairs.
   */

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    auto const spair_range_cp = spair_range;

    /* not necessarily good lower bounds for member variable sizes */
    m_top_rows.reserve(spair_range.size());
    m_bottom_rows.reserve(spair_range.size());
    m_top_coeffs.reserve(spair_range.size());
    m_bottom_coeffs.reserve(spair_range.size());
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);

    std::vector<index_type> gens_ind;

    while (not spair_range.empty())
    {
        spair_monomial_type const lcm = spair_range[0].lcm;

        /* get the range of spair with minimal lcm that are not yet processed;
         * the the spairs are pre-sorted since the range has constant degree */
        auto const lcm_end = std::find_if_not(
            std::cbegin(spair_range) + 1, std::cend(spair_range),
            [lcm](spair_type_t const& sp) { return sp.lcm == lcm; });

        ssize_t const num_lcm =
            std::distance(std::cbegin(spair_range), lcm_end);
        /* each spair contributes with two generators */
        gens_ind.resize(2 * static_cast<size_t>(num_lcm));

        /* copy generators indices with the same lcm into vector */
        std::for_each(std::cbegin(spair_range), lcm_end,
                      [&gens_ind, i = 0ULL](spair_type_t const& sp) mutable {
                          gens_ind[i++] = sp.idx1;
                          gens_ind[i++] = sp.idx2;
                      });

        /* sort all generators with the same lcm and and make them unique
         * (generators found earlier come first) */
        std::sort(std::begin(gens_ind), std::end(gens_ind));

        auto const new_end =
            std::unique(std::begin(gens_ind), std::end(gens_ind));
        gens_ind.erase(new_end, std::cend(gens_ind));

        /* sparsest generator (or smallest by degree) goes in the first
         * position, that is, into the top rows */
        std::ranges::partial_sort(gens_ind, std::begin(gens_ind) + 1, {},
                                  [&basis](index_type const idx) {
#if USE_SPARSEST_REDUCER
                                      return std::make_tuple(basis.length(idx),
                                                             basis.degree(idx));
#else
                                      return std::make_tuple(basis.degree(idx),
                                                             basis.length(idx));
#endif
                                  });

        // TODO thread, this block as well as each iteration of the loop can be
        // safely executed by independent threads
        {
            index_type const idx = gens_ind[0];

            /* the sparsest generator is a reducer and goes to the top part */
            auto const& monomials = basis.monomials(idx);

            m_top_rows.emplace_back(
                create_multiplied_poly_matrix_row(monomials, lcm));

            /* save the generator data to access the coefficients later on */
            coeff_type_ptr const cfs = basis.coefficients(idx).data();

            m_top_coeffs.emplace_back(cfs);
        }

        /* the rest go to the bottom part to be reduced */
        for (size_t i = 1; i < gens_ind.size(); ++i)
        {
            index_type const idx = gens_ind[i];

            auto const& monomials = basis.monomials(idx);

            m_bottom_rows.emplace_back(
                create_multiplied_poly_matrix_row_prefetch(monomials, lcm));

            coeff_type_ptr const cfs = basis.coefficients(idx).data();

            m_bottom_coeffs.emplace_back(cfs);
        }

        /* advance beginning of range & clear used memory */
        spair_range.advance(num_lcm);
        gens_ind.clear();
    }

    /* decrese spair count for each generator in a selected spair before
     * modifying the range */
    std::ranges::for_each(spair_range_cp, [&basis](spair_type_t const& sp) {
        basis.spair_count(sp.idx1)--;
        basis.spair_count(sp.idx2)--;
    });

    /* every time spair_count are decreased new deleted gens. can be created */
    basis.update_deleted_gens();

    /* remove inserted spairs from spair queue */
    spairs.queue.erase(std::begin(spair_range_cp), std::end(spair_range_cp));

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.matrix_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.matrix_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

/*
 * creates the matrix row by multipliying the monomials is poly_mons by the
 * quotient between lm2 and the leading monomial of poly_mons
 */
template <class CoefficientType, class MonomialOrder>
template <class OtherMonomialContext>
typename matrix_f4<CoefficientType, MonomialOrder>::monomial_vect_type
matrix_f4<CoefficientType, MonomialOrder>::create_multiplied_poly_matrix_row(
    basis_monomial_vect_type const& gen_mons,
    monomial<OtherMonomialContext, MonomialOrder> const& mon)
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

template <class CoefficientType, class MonomialOrder>
template <class OtherMonomialContext>
typename matrix_f4<CoefficientType, MonomialOrder>::monomial_vect_type
matrix_f4<CoefficientType, MonomialOrder>::
    create_multiplied_poly_matrix_row_prefetch(
        basis_monomial_vect_type const& gen_mons,
        monomial<OtherMonomialContext, MonomialOrder> const& mon)
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

template <class CoefficientType, class MonomialOrder>
typename matrix_f4<CoefficientType, MonomialOrder>::monomial_vect_type
matrix_f4<CoefficientType, MonomialOrder>::create_poly_matrix_row(
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

template <class CoefficientType, class MonomialOrder>
void matrix_f4<CoefficientType, MonomialOrder>::symbolic_preprocessing(
    polynomial_basis<CoefficientType, MonomialOrder> const& basis)
{
    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

#if USE_DELETED_REDUCERS
    size_t const num_gens = basis.num_gens();

    std::vector<index_type> gens_ind(num_gens);
    std::vector<basis_monomial_type> gens_lm(num_gens);
    aligned_vector<divmask_type> gens_sdm(num_gens);

    /* pack data for reducer search */
    for (size_t gen_idx = 0; gen_idx < num_gens; ++gen_idx)
    {
        gens_ind[gen_idx] = static_cast<index_type>(gen_idx);
        gens_lm[gen_idx]  = basis.lead_mon(gen_idx);
        gens_sdm[gen_idx] = basis.lead_sdm(gen_idx);
    }
#else
    size_t const num_gens = basis.num_nondel_gens();

    std::vector<index_type> gens_ind(num_gens);
    std::vector<basis_monomial_type> gens_lm(num_gens);
    aligned_vector<divmask_type> gens_sdm(num_gens);

    /* pack data for reducer search */
    for (size_t gen_idx = 0, i = 0; gen_idx < basis.num_gens(); ++gen_idx)
    {
        size_t const nondel_idx = basis.nondel_gen_idx(gen_idx);

        /* generator is marked as deleted */
        if (nondel_idx == std::numeric_limits<index_type>::max())
            continue;

        gens_ind[i]   = static_cast<index_type>(gen_idx);
        gens_lm[i]    = basis.lead_mon(gen_idx);
        gens_sdm[i++] = basis.lead_sdm(gen_idx);
    }
#endif

    /* sort generators data by sparsest generator first & break ties with
     * degree; or sort generators by degree and break ties by length */
    std::ranges::stable_sort(std::views::zip(gens_ind, gens_lm, gens_sdm), {},
                             [&basis](auto const gen_data) {
                                 index_type const idx = std::get<0>(gen_data);
#if USE_SPARSEST_REDUCER
                                 return std::make_tuple(basis.length(idx),
                                                        basis.degree(idx));
#else
                                 return std::make_tuple(basis.degree(idx),
                                                        basis.length(idx));
#endif
                             });

    /* iterate directly over all the monomials inserted in the matrix hash
     * table; iterators remain valid after rehashing */
    for (auto it = m_mon_set.cbegin(); it != m_mon_set.cend(); ++it)
    {
        monomial_type const mon = *it;

        /* the monomial is an spair lcm and already has a reducer in the matrix;
           notice that for new reducers the lm was processed previously */
        if (mon.data().idx == 1)
            continue;

        size_t const idx = find_multiplied_reducer(mon, gens_lm, gens_sdm);

        /* no reducer has been found */
        if (idx == num_gens)
            continue;

        index_type const reducer_idx = gens_ind[idx];

        auto const& mons_red = basis.monomials(reducer_idx);
        auto const cfs_red   = basis.coefficients(reducer_idx).data();

        /* create multiplied poly. from reducer and insert it into the matrix */
        auto const& new_mons =
            create_multiplied_poly_matrix_row_prefetch(mons_red, mon);

        m_top_rows.emplace_back(new_mons);
        m_top_coeffs.emplace_back(cfs_red);
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.symbolic_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.symbolic_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    stats.rows_reduced += num_bottom_rows();
}

template <class CoefficientType, class MonomialOrder>
size_t matrix_f4<CoefficientType, MonomialOrder>::find_multiplied_reducer(
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

template <class CoefficientType, class MonomialOrder>
void matrix_f4<CoefficientType, MonomialOrder>::convert_monomials_to_columns(
    gamba_params const& params)
{
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

    /* There is no need to sort the individual rows since the rows/generators
     * are sorted by the given monomial order and the stable partition preserves
     * this order */

    /* sort top rows by pivot order */
    std::ranges::sort(
        std::views::zip(m_top_rows, m_top_coeffs),
        [](index_type const lhs, index_type const rhs) { return lhs < rhs; },
        [](auto const& row_cfs) { return std::get<0>(row_cfs)[0].data().idx; });

    /* sort bottom rows by (reverse) pivot order */
    std::ranges::sort(
        std::views::zip(m_bottom_rows, m_bottom_coeffs),
        [](index_type const lhs, index_type const rhs) { return lhs > rhs; },
        [](auto const& row_cfs) { return std::get<0>(row_cfs)[0].data().idx; });

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.convert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.convert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    if (params.rounds_info)
    {
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

        std::cout << std::format("{:>10} x {:<-9} {:>6.2f}%",
                                 m_top_rows.size() + m_bottom_rows.size(),
                                 m_mon_set.size(), density)
                  << std::flush;
    }
}

template <class CoefficientType, class MonomialOrder>
void matrix_f4<CoefficientType, MonomialOrder>::
    convert_monomials_to_columns_reduce(gamba_params const& params)
{
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

    /* helper lambda */
    auto row_indices = []<typename T>(T& row) {
        using ptr_type = std::conditional_t<std::is_const_v<T>,
                                            index_type const*, index_type*>;

        return reinterpret_cast<ptr_type>(row.data());
    };

    for (auto& row : m_top_rows)
    {
        std::transform(std::cbegin(row), std::cend(row), row_indices(row),
                       [](monomial_type const mon) { return mon.data().idx; });
    }

    for (auto& row : m_bottom_rows)
    {
        std::transform(std::cbegin(row), std::cend(row), row_indices(row),
                       [](monomial_type const mon) { return mon.data().idx; });
    }

    /* there is no need to sort the individual rows since the rows/generators
     * are sorted by the given monomial order and the stable partition preserves
     * this order */

    /* sort top rows by pivot order */
    std::ranges::sort(
        std::views::zip(m_top_rows, m_top_coeffs),
        [](index_type const lhs, index_type const rhs) { return lhs < rhs; },
        [row_indices](auto const& row_cfs) {
            return row_indices(std::get<0>(row_cfs))[0];
        });

    /* sort bottom rows by pivot order */
    std::ranges::sort(
        std::views::zip(m_bottom_rows, m_bottom_coeffs),
        [](index_type const lhs, index_type const rhs) { return lhs > rhs; },
        [row_indices](auto const& row_cfs) {
            return row_indices(std::get<0>(row_cfs))[0];
        });

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.convert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.convert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    if (params.rounds_info)
    {
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

        std::cout << std::format("{:>10} x {:<-9} {:>6.2f}%",
                                 m_top_rows.size() + m_bottom_rows.size(),
                                 m_mon_set.size(), density)
                  << std::flush;
    }
}

template <class CoefficientType, class MonomialOrder>
void matrix_f4<CoefficientType, MonomialOrder>::insert_generators_reduce(
    polynomial_basis<CoefficientType, MonomialOrder> const& basis)
{
    size_t const num_reduced_gens = basis.reduced_indices().size();

    m_bottom_rows.reserve(num_reduced_gens);
    m_bottom_coeffs.reserve(num_reduced_gens);
    m_new_coefs.reserve(num_reduced_gens);
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);

    for (size_t const i : basis.reduced_indices())
    {
        auto const& monomials = basis.monomials(i);
        m_bottom_rows.emplace_back(create_poly_matrix_row(monomials));

        /* if using reduced generators the leading monomial is already done
         * for symbolic preprocessing */
        monomial_type const lm_row = m_bottom_rows.back()[0];
        lm_row.data().idx          = 1;

        coeff_type_ptr const cfs = basis.coefficients(i).data();
        m_bottom_coeffs.emplace_back(cfs);
    }
}

template <class CoefficientType, class MonomialOrder>
void matrix_f4<CoefficientType, MonomialOrder>::extract_new_rows(
    linalg_v3<CoefficientType, MonomialOrder>& linalg)
{
    using row_type = linalg_v3<CoefficientType, MonomialOrder>::row_type;

    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    // m_top_rows.clear();
    // m_top_rows.shrink_to_fit();
    // m_bottom_rows.clear();
    // m_bottom_rows.shrink_to_fit();

    /* reorder the matrix cols. -> monomials map so that the relative non-pivot
     * column indices can be mapped directly to the corresponding monomials */
    std::stable_partition(std::begin(m_col_to_mon), std::end(m_col_to_mon),
                          [&linalg](monomial_type const mon) {
                              return not linalg.m_is_pivot[mon.data().idx];
                          });

    for (size_t i = 0; i < linalg.m_newrows_buffers.size(); ++i)
    {
        auto const& new_rows = linalg.m_newrows_buffers[i];
        auto const& new_pivs = linalg.m_newrows_piv_ind[i];

        for (size_t row_idx = 0; row_idx < new_pivs.size(); ++row_idx)
        {
            index_type const col_idx = new_pivs[row_idx];

            assert(row_idx < NUM_ROWS_BUFFER);
            assert(col_idx < linalg.m_num_nonpiv_cols);

            /* add 1 to num_nnz since the pivot term is not counted */
            size_t const abs_idx = NUM_ROWS_BUFFER * i + row_idx;
            size_t const num_nnz = linalg.m_newrows_num_nnz[abs_idx] + 1;

            auto const [new_mon, new_cfs] =
                basis_type::allocate_polynomial(num_nnz);

            /* use monomial array to store matrix column indices */
            auto* const new_ind = reinterpret_cast<index_type*>(new_mon);

            m_new_rows.emplace_back(new_ind, num_nnz);
            m_new_coefs.emplace_back(new_cfs, num_nnz);

            /* new rows are already monic so first term is equal to 1 * r */
            new_ind[0] = col_idx;
            new_cfs[0] = m_field.r;

            length_type new_size = 1;

            row_type const* const row_ptr = new_rows.data() + row_idx;

            for (size_t j = col_idx + 1; j < linalg.m_num_nonpiv_cols; ++j)
            {
                assert(row_ptr[j * NUM_ROWS_BUFFER]
                       < std::numeric_limits<coeff_type>::max());

                row_type const c = row_ptr[j * NUM_ROWS_BUFFER];

                /* coefficients in new rows buffers are already normalized */
                if (c == 0)
                    continue;

                new_ind[new_size] = static_cast<index_type>(j);
                new_cfs[new_size] = static_cast<coeff_type>(c);
                new_size++;
            }

            assert(linalg.m_num_bottom_rows >= m_new_coefs.size());

            assert(new_cfs[0] == m_field.r);
            assert(new_cfs[new_size - 1] != 0);
        }

        /* free memory of current buffer before allocating memory for next
         * buffer's rows; hopefully this helps the allocator reuse memory */
        linalg.m_newrows_buffers[i].clear();
        linalg.m_newrows_buffers[i].shrink_to_fit();
        delete[] linalg.m_newrows_piv_ind[i].data();
    }

    linalg.m_newrows_buffers.clear();
    linalg.m_newrows_piv_ind.clear();
    linalg.m_newrows_is_pivot.clear();
    linalg.m_newrows_num_nnz.clear();

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.insert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.insert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

template <class CoefficientType, class MonomialOrder>
double matrix_f4<CoefficientType, MonomialOrder>::memory_usage() const
{
    double mem_size = 0.0;

    mem_size += memory_size(m_top_rows);
    for (auto const& row : m_top_rows)
        mem_size += memory_size(row);

    mem_size += memory_size(m_bottom_rows);
    for (auto const& row : m_bottom_rows)
        mem_size += memory_size(row);

    mem_size += memory_size(m_top_coeffs);
    mem_size += memory_size(m_bottom_coeffs);

    mem_size += memory_size(m_col_to_mon);

    mem_size += m_mon_set.memory_usage();

    mem_size += memory_size(m_new_coefs);
    for (auto const& cfs : m_new_coefs)
        mem_size += memory_size(cfs);

    mem_size += memory_size(m_new_rows);
    for (auto const& row : m_new_rows)
        mem_size += memory_size(row);

    return mem_size;
}

template <class CoefficientType, class MonomialOrder>
void matrix_f4<CoefficientType, MonomialOrder>::clear()
{
    m_top_rows.clear();
    m_bottom_rows.clear();

    m_top_coeffs.clear();
    m_bottom_coeffs.clear();

    m_mon_set.clear();

    m_col_to_mon.clear();

    m_new_rows.clear();
    m_new_coefs.clear();
}

}  // namespace gamba
