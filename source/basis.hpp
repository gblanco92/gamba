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
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "config.hpp"
#include "container.hpp"
#include "divmask.hpp"
#include "field.hpp"
#include "io.hpp"
#include "monomial.hpp"
#include "order.hpp"
#include "stats.hpp"

namespace gamba
{

template <class CoefficientType, class MonomialOrder>
class polynomial_basis
{
public:
    using coeff_type     = CoefficientType;
    using coeff_type_ptr = coeff_type const*;
    using monomial_order = MonomialOrder;

    using monomial_context = basis_hashtable;
    using monomial_type    = monomial<monomial_context, monomial_order>;
    using monomial_init    = typename monomial_type::monomial_init_t;

    using hash_type      = typename monomial_type::hash_type;
    using degree_type    = typename monomial_type::degree_type;
    using index_type     = typename monomial_type::index_type;
    using index_type_ptr = index_type const*;

    using length_type = uint32_t;
    using count_type  = uint32_t;

    using divmap_type = divmask_map<monomial_order, monomial_type>;

    using matrix_monomial_context = matrix_hashtable;
    using matrix_monomial_type =
        monomial<matrix_monomial_context, monomial_order>;

    using monomial_vect_type = std::span<monomial_type>;
    using index_vect_type    = std::span<index_type>;
    using coeff_vect_type    = std::span<coeff_type>;

private:
    static constexpr auto is_divisible =
        monomial_type::template check_monomial_division_divmask<
            monomial_context>;

public:
    polynomial_basis(size_t const n_vars, uint32_t const fld_chr);

    ~polynomial_basis();

    /* make this class non-copyable */
    polynomial_basis(polynomial_basis const&)            = delete;
    polynomial_basis& operator=(polynomial_basis const&) = delete;

    void clear();

    /* access member functions */
    size_t num_gens() const { return m_num_gens; }

    length_type length(size_t const i) const
    {
        return static_cast<length_type>(m_mons[i].size());
    }

    degree_type degree(size_t const i) const { return m_degs[i]; }

    monomial_type lead_mon(size_t const i) const { return m_mons[i][0]; }

    divmask_type lead_sdm(size_t const i) const { return m_lead_sdm[i]; }

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

    bool is_trivial() const { return m_trivial; }

    divmap_type const& divmap() const { return m_divmap; }

    monomial_vect_type const& monomials(size_t const i) const
    {
        return m_mons[i];
    }

    coeff_vect_type const& coefficients(size_t const i) const
    {
        return m_coefs[i];
    }

    size_t divmasks_version() const { return m_divmasks_version; }

    count_type& spair_count(size_t const i) { return m_spair_count[i]; }

    static std::pair<monomial_type*, coeff_type*> allocate_polynomial(
        size_t const num_terms);

    void import_generators(generators_data const& data);

    /* after an update process enforce the following invariant:
     *  redundant[reduced_lm[i]] == false */
    void update_reduced_gens(size_t const prev_num_gens);

    template <class MatrixType>
    void insert_new_rows_reduce(MatrixType& matrix);

    template <class MatrixType2>
    void insert_new_rows_echelon(MatrixType2& matrix);

    void remove_redundant_gens();

    void update_deleted_gens();

    /* I/O and info member functions */
#if defined DEBUG
    void write_generators(std::ostream& os) const;
#endif

    void export_generators(generators_data& data) /*const*/;

    void print_info(std::ostream& os) const;

    double memory_usage() const;

private:
    void import_generator(size_t const i, generators_data const& data);

    void import_matrix_row(index_vect_type& row,
                           coeff_vect_type const& cfs,
                           std::vector<matrix_monomial_type> const& col_to_mon);

public:
    /* number of variables in the polynomial ring */
    size_t const num_vars;

    /* field metadata */
    field_traits<coeff_type> const field;

private:
    /* numer of total generators (non-reduced) in the basis */
    size_t m_num_gens{};

    /* index of the hashed monomials in the monomial set */
    std::vector<monomial_vect_type> m_mons{};

    /* coefs[i][j] corresponds to the coeff. of the mon. index by mons[i][j] */
    std::vector<coeff_vect_type> m_coefs{};

    /* unordered (flat) set storing all monomials appearing in the basis */
    monomial_set<monomial_type> m_mon_set{};

    /* total degree of each polynomial in the basis */
    std::vector<degree_type> m_degs{};

    /* divisibility mask for each generator's leading monomial */
    std::vector<divmask_type> m_lead_sdm{};

    /* helper class to generate the divisibility masks */
    divmap_type m_divmap;

    /* version of divmap used in the computation of current lead_sdm */
    size_t m_divmasks_version{0UL};

    /* polynomials of the basis made redundant by G-M update; make it aligned so
     * it can be updated concurrently, avoid std::vector<bool> specialization */
    aligned_vector<uint8_t> m_redundant{};

    /* keeps the index of top reduced polynomials in the basis; this means that
     * for all index i: redundant[reduced_lm[i]] == false for all i */
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

    /* the ideal is trivial */
    bool m_trivial{false};
};

template <class CoefficientType, class MonomialOrder>
polynomial_basis<CoefficientType, MonomialOrder>::polynomial_basis(
    size_t const n_vars,
    uint32_t const fld_chr) :
        num_vars{n_vars},

        field{fld_chr},

        m_divmap{n_vars}
{}

template <class CoefficientType, class MonomialOrder>
polynomial_basis<CoefficientType, MonomialOrder>::~polynomial_basis()
{
    for (auto& gen_mons : m_mons)
        delete[] gen_mons.data();
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::clear()
{
    m_num_gens = 0UL;

    for (auto& gen_mons : m_mons)
        delete[] gen_mons.data();

    m_mons.clear();
    m_mons.shrink_to_fit();

    m_coefs.clear();
    m_coefs.shrink_to_fit();

    m_mon_set.clear();

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

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::import_generators(
    generators_data const& data)
{ /*
   * The import process consists of the following steps:
   * 1. Allocate memory for all the generators.
   * 2. Insert monomials to hash tables & store corresponding coefficient.
   * 3. Sort monomials & coefficients with respect to the monomial order.
   * 4. Sum and remove equal terms and monomials.
   * 5. Normalize the coefficients of each polynomial.
   * 6. (Char = p) Transforms coefficients to modular representation.
   * 7. Compute degree of each polynomial generator.
   * 8. Remove empty, that is equal to zero, generators.
   * 9. Sort generators by leading monomial.
   * 10. Check if the genereting set is homogeneous.
   * 11. Generate divisibility masks for all monomials in the basis.
   */
    m_num_gens = data.num_gens;

    m_coefs.resize(m_num_gens);
    m_mons.resize(m_num_gens);
    m_degs.resize(m_num_gens);
    m_mon_set.reserve(HASHTABLE_INIT_SIZE);

    for (size_t idx = 0; idx < m_num_gens; ++idx)
        import_generator(idx, data);

    /* remove empty generators */
    auto generators           = std::views::zip(m_mons, m_coefs, m_degs);
    auto const [new_end, end] = std::ranges::remove(
        generators, 0, [](auto const& p) { return get_mons(p).size(); });

    auto const gens_begin = std::ranges::cbegin(generators);

    /* deallocate before removal */
    for (auto const& gen : std::ranges::subrange{new_end, end})
        delete[] get_mon(gen).data();

    /* do the actual removal */
    m_mons.erase(std::cbegin(m_mons) + (new_end - gens_begin),
                 std::cend(m_mons));
    m_coefs.erase(std::cbegin(m_coefs) + (new_end - gens_begin),
                  std::cend(m_coefs));
    m_degs.erase(std::cbegin(m_degs) + (new_end - gens_begin),
                 std::cend(m_degs));

    /* true number of generators after removing empty generators */
    m_num_gens = m_mons.size();

    /* resize directly data that is constant for all new generators */
    m_redundant.resize(m_num_gens, false);
    m_spair_count.resize(m_num_gens, 0);

    m_num_nondel_gens = m_num_gens;

    /* in the begining no generator is marked as deleted */
    m_nondel_gens.resize(m_num_gens);
    std::iota(std::begin(m_nondel_gens), std::end(m_nondel_gens), 0U);

    generators = std::views::zip(m_mons, m_coefs, m_degs);

    /* sort generators by leading monomial; this ordering will change how the
     * Gebauer-Moeller process creates and updates spairs */
    std::ranges::sort(generators, [](auto const& x, auto const& y) {
#if INSERT_ELEMENTS_DECREASING == 1
        return monomial_order{}(get_mons(x)[0], get_mons(y)[0]) > 0;
#else
        return monomial_order{}(get_mons(x)[0], get_mons(y)[0]) < 0;
#endif
    });

    /* helper lambda to check if a single polynomials is homogeneous */
    auto const is_homogeneous = [](auto const& poly) {
        return std::ranges::all_of(get_mons(poly),
                                   [&poly](monomial_type const mon) {
                                       return mon.degree() == get_deg(poly);
                                   });
    };

    /* check if the basis is homogeneous after all the cleaning */
    m_is_homogeneous = std::ranges::all_of(generators, is_homogeneous);

    /* generate divisibility masks after importing generators into basis */
    m_divmap.update_map(m_mon_set);

    m_lead_sdm.resize(m_num_gens);
    /* compute divmask for each leading term in new basis */
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        monomial_type const lm = m_mons[i][0];
        m_lead_sdm[i]          = m_divmap.compute_divmask(lm);
    }

    m_divmasks_version = m_divmap.version();
}

template <class CoefficientType, class MonomialOrder>
std::pair<
    typename polynomial_basis<CoefficientType, MonomialOrder>::monomial_type*,
    typename polynomial_basis<CoefficientType, MonomialOrder>::coeff_type*>
polynomial_basis<CoefficientType, MonomialOrder>::allocate_polynomial(
    size_t const num_terms)
{
    constexpr size_t const std_align = __STDCPP_DEFAULT_NEW_ALIGNMENT__;

    size_t const num_bytes_mons  = num_terms * sizeof(monomial_type);
    size_t const num_bytes_coefs = num_terms * sizeof(coeff_type);
    size_t const num_bytes_padd  = padding<std_align>(num_bytes_mons);

    size_t const num_bytes = num_bytes_mons + num_bytes_coefs + num_bytes_padd;
    size_t const offset_bytes = num_bytes_mons + num_bytes_padd;

    /* allocate monomials & coefficients in a single allocation */
    auto* mem_ptr = new char[num_bytes];
    auto* mon_ptr = reinterpret_cast<monomial_type*>(mem_ptr);
    auto* cfs_ptr = reinterpret_cast<coeff_type*>(mem_ptr + offset_bytes);

    /* both memory blocks have the same alignment as if allocated separately */
    assert(reinterpret_cast<std::uintptr_t>(mem_ptr) % std_align == 0);
    assert(reinterpret_cast<std::uintptr_t>(cfs_ptr) % std_align == 0);

    return std::make_pair(mon_ptr, cfs_ptr);
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::import_generator(
    size_t const i,
    generators_data const& data)
{
    size_t const num_terms = data.lens[i];

    auto const [mon_ptr, cfs_ptr] = allocate_polynomial(num_terms);

    m_mons[i]  = {mon_ptr, num_terms};
    m_coefs[i] = {cfs_ptr, num_terms};

    // NOLINTNEXTLINE
    size_t const offset = std::accumulate(&data.lens[0], &data.lens[i], 0ULL);

    for (size_t j = 0; j < num_terms; ++j)
    {
        /* convert input monomial into own format & compute hash value */
        auto const mon = monomial_init::construct(
            &data.exps[num_vars * (offset + j)], num_vars);

        /* insert monomial in basis hash table */
        auto const [it, _] = m_mon_set.insert(mon);

        /* copy monomial into generator vector */
        m_mons[i][j] = *it;
    }

    /* import generators coefficients */
    if constexpr (std::is_same_v<coeff_type, mpz_class>)
    {
        std::copy(&data.coeffs[offset], &data.coeffs[offset + data.lens[i]],
                  m_coefs[i].data());
    }
    else
    {
        std::copy(&data.coeffs_modp[offset],
                  &data.coeffs_modp[offset + data.lens[i]], m_coefs[i].data());
    }

    auto const poly = std::views::zip(m_mons[i], m_coefs[i]);
    /* sort monomials and coefficients with the given (decreasing) mon. order */
    std::ranges::sort(poly, [](auto const x, auto const y) {
        return monomial_order{}(get_mon(x), get_mon(y)) > 0;
    });

    auto const rpoly = poly | std::views::reverse;
    /* add coeffs with same monomials in the first coefficient of the run */
    std::ranges::transform(
        rpoly | std::views::take(poly.size() - 1), rpoly | std::views::drop(1),
        std::rbegin(m_coefs[i]) + 1, [this](auto const x, auto const y) {
            return get_mon(x) == get_mon(y)
                     ? field.add(get_coef(x), get_coef(y))
                     : get_coef(y);
        });

    auto const poly_begin = std::ranges::begin(poly);

    /* remove repeated coefficients/monomials */
    auto const [new_end0, _0] = std::ranges::unique(poly, {}, get_mon);

    /* remove zero coefficients */
    auto const [new_end, _] =
        std::ranges::remove(poly_begin, new_end0, 0, get_coef);

    /* new size after all simplifications */
    auto const new_size =
        static_cast<size_t>(std::distance(poly_begin, new_end));

    /* do the actual removal by simply updating the size */
    m_mons[i]  = {m_mons[i].data(), new_size};
    m_coefs[i] = {m_coefs[i].data(), new_size};

    /* early return if zero polynomial is found */
    if (new_size == 0)
        return;

    /* normalize coefs: if char > 0 make monic, if char = 0 remove content */
    if constexpr (std::is_same_v<coeff_type, mpz_class>)
    {
        mpz_class const gcd0 = std::accumulate(
            std::begin(m_coefs[i]), std::end(m_coefs[i]), mpz_class{0},
            [](mpz_class const& acc, mpz_class const& c) {
                return gcd(acc, c);
            });

        /* remove content from rational coefficients */
        std::ranges::for_each(m_coefs[i], [&gcd0](mpz_class& c) { c /= gcd0; });
    }
    else
    {
        coeff_type const inv = field.inverse(m_coefs[i][0]);

        /* make polynomial monic */
        std::ranges::for_each(m_coefs[i], [inv, this](coeff_type& c) {
            c = field.multiply(c, inv);
        });

        assert(m_coefs[i][0] == 1);
    }

    /* transform mod p coefs to modular space *after* normalization */
    if constexpr (not std::is_same_v<coeff_type, mpz_class>)
    {
        std::transform(std::begin(m_coefs[i]), std::end(m_coefs[i]),
                       std::begin(m_coefs[i]), [this](coeff_type const x) {
                           return field.transform(x);
                       });
    }
    assert(m_coefs[i][0] == field.r);

    /* compute degree after all simplifications */
    if constexpr (is_degree_order_v<monomial_order>)
    {
        m_degs[i] = m_mons[i][0].degree();
    }
    else
    {
        auto const it = std::ranges::max_element(
            m_mons[i], {},
            [](monomial_type const mon) { return mon.degree(); });

        m_degs[i] = it->degree();
    }

    if (m_degs[i] == 0)
        m_trivial = true;
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::update_reduced_gens(
    size_t const prev_num_gens)
{
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

    stats.redundant_elem =
        static_cast<ssize_t>(m_num_gens - m_reduced_gens.size());
}

template <class CoefficientType, class MonomialOrder>
template <class MatrixType>
void polynomial_basis<CoefficientType, MonomialOrder>::insert_new_rows_reduce(
    MatrixType& matrix)
{
    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    assert(m_num_gens == 0);
    assert(m_mons.size() == 0);
    assert(m_coefs.size() == 0);
    assert(m_redundant.empty());

    m_num_gens = matrix.m_bottom_rows.size();

    /* all generators are non-reduced after reduce phase */
    m_redundant.resize(m_num_gens, false);
    m_reduced_gens.resize(m_num_gens);
    std::iota(std::begin(m_reduced_gens), std::end(m_reduced_gens), 0UL);

    /* reserve memory for the new generators */
    m_mons.reserve(m_num_gens);
    m_coefs.reserve(m_num_gens);
    m_degs.reserve(m_num_gens);

    /* rows come presorted from the reduce phase */
    for (size_t i = 0; i < matrix.m_bottom_rows.size(); ++i)
    {
        index_vect_type& row = matrix.m_new_rows[i];
        coeff_vect_type& cfs = matrix.m_new_coefs[i];

        import_matrix_row(row, cfs, matrix.m_col_to_mon);
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.insert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.insert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

template <class CoefficientType, class MonomialOrder>
template <class MatrixType>
void polynomial_basis<CoefficientType, MonomialOrder>::insert_new_rows_echelon(
    MatrixType& matrix)
{
    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    size_t const num_new_rows  = matrix.m_new_rows.size();
    size_t const prev_num_gens = m_num_gens;

    m_num_gens += num_new_rows;

    /* resize directly data that is constant for all new generators */
    m_redundant.resize(m_num_gens, false);
    m_spair_count.resize(m_num_gens, 0);

    /* all newly inserted generators are non-deleted */
    m_nondel_gens.resize(m_num_gens);
    std::iota(std::begin(m_nondel_gens) + static_cast<ssize_t>(prev_num_gens),
              std::end(m_nondel_gens), m_num_nondel_gens);

    m_num_nondel_gens += num_new_rows;

    /* reserve memory for the new generators */
    m_mons.reserve(m_num_gens);
    m_coefs.reserve(m_num_gens);
    m_degs.reserve(m_num_gens);

    /* rows come sorted from the reduce phase by decreasing order of leading
     * monomials (highest leading monomials go first); how these elements are
     * inserted into the basis will change how Gebauer-Moeller updates spairs */
#if INSERT_ELEMENTS_DECREASING == 1
    for (size_t i : std::views::iota(0ULL, num_new_rows))
#else
    for (size_t i : std::views::iota(0ULL, num_new_rows) | std::views::reverse)
#endif
    {
        index_vect_type& row = matrix.m_new_rows[i];
        coeff_vect_type& cfs = matrix.m_new_coefs[i];

        import_matrix_row(row, cfs, matrix.m_col_to_mon);
    }

    matrix.clear();

    /* update div. masks after importing new generators into basis */
    m_divmap.update_map(m_mon_set);

    m_lead_sdm.resize(m_num_gens);

    /* only update sdm of previous generators if new divmap version */
    size_t const start_idx =
        m_divmap.version() != m_divmasks_version ? 0UL : prev_num_gens;

    /* update divmask of leading monomials for all generators since spairs
     * can depend on non-reduced generators */
    for (size_t i = start_idx; i < m_num_gens; ++i)
    {
        monomial_type const lm = m_mons[i][0];
        m_lead_sdm[i]          = m_divmap.compute_divmask(lm);
    }

    m_divmasks_version = m_divmap.version();

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats.insert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats.insert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::import_matrix_row(
    polynomial_basis<CoefficientType, MonomialOrder>::index_vect_type& row,
    polynomial_basis<CoefficientType, MonomialOrder>::coeff_vect_type const&
        cfs,
    std::vector<matrix_monomial_type> const& col_to_mon)
{
    size_t const num_terms = row.size();

    auto const ind_ptr = reinterpret_cast<index_type_ptr>(row.data());
    auto const mon_ptr = reinterpret_cast<monomial_type*>(row.data());

    monomial_vect_type mons = {mon_ptr, num_terms};

    for (size_t i = 0; i < num_terms; ++i)
    {
        /* undo the columns <-> monomials transformation */
        matrix_monomial_type mat_mon = col_to_mon[ind_ptr[i]];

        /* copy matrix monomial into a basis monomial */
        auto const mon = monomial_init::copy(mat_mon);

        /* insert basis monomial into basis hash table */
        auto const [it, _] = m_mon_set.insert(mon);

        /* copy monomial to generator vector */
        mons[i] = *it;
    }

    m_mons.emplace_back(mons);
    m_coefs.emplace_back(cfs);

    /* compute degree of new generator */
    if constexpr (is_degree_order_v<monomial_order>)
    {
        m_degs.emplace_back(m_mons.back()[0].degree());
    }
    else
    {
        auto const it = std::ranges::max_element(
            m_mons.back(), {},
            [](monomial_type const mon) { return mon.degree(); });

        m_degs.emplace_back(it->degree());
    }

    if (m_degs.back() == 0)
        m_trivial = true;
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::remove_redundant_gens()
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

    stats.redundant_elem =
        static_cast<ssize_t>(m_num_gens - m_reduced_gens.size());
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::update_deleted_gens()
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

#if defined DEBUG
template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::write_generators(
    std::ostream& os) const
{
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        os << "#" << i << ": deg = " << m_degs[i] << ", ";

        for (size_t j = 0; j < m_mons[i].size(); ++j)
        {
            std::string const& mon_str = monomial2string(m_mons[i][j]);

            os << (j != 0 ? " + " : "") << static_cast<uint64_t>(m_coefs[i][j]);
            os << (mon_str.empty() ? "" : "*") << mon_str;
        }

        os << std::endl;
    }
}
#endif

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::export_generators(
    generators_data& data) /*const*/
{
    using var_type     = generators_data::var_type;
    using coeff_p_type = generators_data::coeff_p_type;

    size_t const num_red_gens = m_reduced_gens.size();

    data.num_vars   = static_cast<var_type>(num_vars);
    data.field_char = static_cast<var_type>(field.n);
    data.num_gens   = static_cast<var_type>(num_red_gens);

    size_t const num_mons = std::accumulate(
        std::cbegin(m_mons), std::cend(m_mons), 0ULL,
        [](size_t acc, auto const& v) { return acc + v.size(); });

    data.coeffs_modp.resize(num_mons);
    // data.coeffs.resize(num_mons); // TODO: char = 0
    data.exps.resize(num_mons * num_vars);
    data.lens.resize(data.num_gens);

    /* write generators in increasing order of leading monomial */
    for (size_t i = 0, offset = 0; i < num_red_gens; ++i)
    {
        size_t const idx = m_reduced_gens[i];

        /* transform mod p coefs from modular space */
        if constexpr (not std::is_same_v<coeff_type, mpz_class>)
        {
            std::transform(std::begin(m_coefs[idx]), std::end(m_coefs[idx]),
                           std::begin(m_coefs[idx]),
                           [this](coeff_type const x) {
                               return field.reduce_normalize(x);
                           });
        }

        for (size_t j = 0; j < m_mons[idx].size(); ++j)
        {
            monomial2exponent(m_mons[idx][j],
                              &data.exps[(offset + j) * num_vars]);

            data.coeffs_modp[offset + j] =
                static_cast<coeff_p_type>(m_coefs[idx][j]);
        }

        data.lens[i] = m_mons[idx].size();

        offset += m_mons[idx].size();
    }
}

template <class CoefficientType, class MonomialOrder>
void polynomial_basis<CoefficientType, MonomialOrder>::print_info(
    std::ostream& os) const
{
    size_t const num_mons = std::accumulate(
        std::cbegin(m_mons), std::cend(m_mons), 0ULL,
        [](size_t acc, auto const& v) { return acc + v.size(); });

    auto const max_degree_poly =
        std::max_element(std::cbegin(m_degs), std::cend(m_degs));

    size_t const max_degree = (max_degree_poly != std::cend(m_degs)
                                   ? *max_degree_poly
                                   : std::numeric_limits<degree_type>::max());

    auto const dnum_gens = static_cast<double>(m_num_gens);
    auto const dnum_mon  = static_cast<double>(num_mons);

    os << "************ BASIS INFO ************" << std::endl;
    os << std::setw(25) << std::left << "num. variables:" << num_vars
       << std::endl;
    os << std::setw(25) << std::left
       << "field charac.:" << static_cast<uint64_t>(field.n) << std::endl;
    os << std::setw(25) << std::left << "num. generators:" << m_num_gens
       << std::endl;
    os << std::setw(25) << std::left << "num. monomials:" << num_mons
       << std::endl;
    os << std::setw(25) << std::left << "avg. length:" << dnum_mon / dnum_gens
       << std::endl;
    os << std::setw(25) << std::left << "max. degree:" << max_degree
       << std::endl;
    os << std::setw(25) << std::left << "homogeneous:" << std::boolalpha
       << m_is_homogeneous << std::endl;
    os << "************************************" << std::endl;
}

template <class CoefficientType, class MonomialOrder>
double polynomial_basis<CoefficientType, MonomialOrder>::memory_usage() const
{
    double mem_size = 0.0;

    mem_size += memory_size(m_coefs);
    mem_size += memory_size(m_mons);
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        mem_size += memory_size(m_coefs[i]);
        mem_size += memory_size(m_mons[i]);
    }

    mem_size += m_mon_set.memory_usage();

    mem_size += memory_size(m_degs);

    mem_size += memory_size(m_redundant);

    mem_size += memory_size(m_reduced_gens);

    mem_size += memory_size(m_spair_count);

    mem_size += memory_size(m_nondel_gens);

    return mem_size;
}

}  // namespace gamba
