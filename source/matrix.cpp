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

#include "matrix.hpp"

namespace gamba
{

void matrix_f4::insert_spairs(spair_set& spairs,
                              size_t const num_spairs,
                              base_polynomial_basis& basis)
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

    /* not necessarily good lower bounds for member variable sizes */
    m_top_rows.reserve(num_spairs);
    m_bottom_rows.reserve(num_spairs);
    m_top_coefs.reserve(num_spairs);
    m_bottom_coefs.reserve(num_spairs);
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);

    auto spairs_begin = std::cbegin(spairs.queue);
    auto const spairs_end =
        std::cbegin(spairs.queue) + static_cast<ssize_t>(num_spairs);

    std::vector<index_type> gens_ind;

    while (spairs_begin != spairs_end)
    {
        spair_monomial_type const lcm = spairs_begin->lcm;

        /* get the range of spair with minimal lcm that are not yet processed;
         * the the spairs are pre-sorted since the range has constant degree */
        auto const lcm_end = std::find_if_not(
            spairs_begin + 1, spairs_end,
            [lcm](spair_type const& sp) { return sp.lcm == lcm; });

        ssize_t const num_lcm = std::distance(spairs_begin, lcm_end);
        /* each spair contributes with two generators */
        gens_ind.resize(2 * static_cast<size_t>(num_lcm));

        /* copy generators indices with the same lcm into vector */
        std::for_each(spairs_begin, lcm_end,
                      [&gens_ind, i = 0ULL](spair_type const& sp) mutable {
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
            basis_monomial_vect_type const& monomials = basis.monomials(idx);

            m_top_rows.emplace_back(
                create_multiplied_poly_matrix_row_prefetch(monomials, lcm));

            /* save the generator data to access the coefficients later on */
            void_ptr_type const cfs = basis.v_coefficients(idx);

            m_top_coefs.emplace_back(cfs);
        }

        /* the rest go to the bottom part to be reduced */
        for (size_t i = 1; i < gens_ind.size(); ++i)
        {
            index_type const idx = gens_ind[i];

            basis_monomial_vect_type const& monomials = basis.monomials(idx);

            m_bottom_rows.emplace_back(
                create_multiplied_poly_matrix_row_prefetch(monomials, lcm));

            void_ptr_type const cfs = basis.v_coefficients(idx);

            m_bottom_coefs.emplace_back(cfs);
        }

        /* advance iterator */
        spairs_begin += num_lcm;
        /* clear used memory */
        gens_ind.clear();
    }

    /* decrese spair count for each generator in a selected spair before
     * modifying the range */
    std::for_each(std::cbegin(spairs.queue), spairs_end,
                  [&basis](spair_type const& sp) {
                      basis.spair_count(sp.idx1)--;
                      basis.spair_count(sp.idx2)--;
                  });

    /* every time spair_count are decreased new deleted gens. can be created */
    basis.update_deleted_gens();

    /* remove inserted spairs from spair queue */
    spairs.queue.erase(std::begin(spairs.queue), spairs_end);

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::matrix_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::matrix_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

void matrix_f4::symbolic_preprocessing(base_polynomial_basis const& basis)
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

        basis_monomial_vect_type const& mons_red = basis.monomials(reducer_idx);
        void_ptr_type const cfs_red = basis.v_coefficients(reducer_idx);

        /* create multiplied poly. from reducer and insert it into the matrix */
        monomial_vect_type const& new_mons =
            create_multiplied_poly_matrix_row_prefetch(mons_red, mon);

        m_top_rows.emplace_back(new_mons);
        m_top_coefs.emplace_back(cfs_red);
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::symbolic_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::symbolic_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;

    stats::rows_reduced += num_bottom_rows();
}

void matrix_f4::insert_generators_reduce(base_polynomial_basis const& basis)
{
    size_t const num_reduced_gens = basis.reduced_indices().size();

    m_bottom_rows.reserve(num_reduced_gens);
    m_bottom_coefs.reserve(num_reduced_gens);
    m_new_coefs.reserve(num_reduced_gens);
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);

    for (size_t const i : basis.reduced_indices())
    {
        basis_monomial_vect_type const& monomials = basis.monomials(i);
        m_bottom_rows.emplace_back(create_poly_matrix_row(monomials));

        /* if using reduced generators the leading monomial is already done
         * for symbolic preprocessing */
        monomial_type const lm_row = m_bottom_rows.back()[0];
        lm_row.data().idx          = 1;

        void_ptr_type const cfs = basis.v_coefficients(i);

        m_bottom_coefs.emplace_back(cfs);
    }
}

double matrix_f4::memory_usage() const
{
    double mem_size = 0.0;

    mem_size += memory_size(m_top_rows);
    for (auto const& row : m_top_rows)
        mem_size += memory_size(row);

    mem_size += memory_size(m_bottom_rows);
    for (auto const& row : m_bottom_rows)
        mem_size += memory_size(row);

    mem_size += memory_size(m_top_coefs);
    mem_size += memory_size(m_bottom_coefs);

    mem_size += memory_size(m_col_to_mon);

    mem_size += m_mon_set.memory_usage();

    mem_size += memory_size(m_new_coefs);
    mem_size += memory_size(m_new_rows);

    /* approx. coeff_type by uint32_t */
    for (auto const& row : m_new_rows)
        mem_size += 2 * memory_size(row);

    return mem_size;
}

void matrix_f4::clear()
{
    m_top_rows.clear();
    m_bottom_rows.clear();

    m_top_coefs.clear();
    m_bottom_coefs.clear();

    m_mon_set.clear();

    m_col_to_mon.clear();

    m_new_rows.clear();
    m_new_coefs.clear();
}

}  // namespace gamba
