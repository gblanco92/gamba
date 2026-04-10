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
#include <memory>
#include <ranges>
#include <vector>

#include "alloc_poly.hpp"
#include "base_basis.hpp"
#include "config.hpp"
#include "field.hpp"
#include "io.hpp"
#include "learn.hpp"
#include "logger.hpp"
#include "matrix.hpp"
#include "order.hpp"
#include "stats.hpp"
#include "thirdparty/flintxx.hpp"
#include "utils.hpp"

namespace gamba
{

template <class CoefficientType>
class polynomial_basis : public base_polynomial_basis
{
public:
    using coeff_type      = CoefficientType;
    using coeff_ptr_type  = coeff_type const*;
    using coeff_vect_type = std::span<coeff_type>;

    polynomial_basis(size_t const n_vars, uint32_t const fld_chr) :
            base_polynomial_basis{n_vars, fld_chr},
            field{fld_chr}
    {}

    /* make all template specialization friends of each other so we can access
     * the private construtor below */
    template <class OtherCoefficientType>
    friend class polynomial_basis;

private:
    /* this constructor *must* be private since the resulting basis is partially
     * constructed; only used when reducing a QQ basis modulo a prime number */
    explicit polynomial_basis(uint32_t const prime) :
            base_polynomial_basis(prime),  // only sets m_field_char = prime
            field{prime}
    {}

public:
    ~polynomial_basis() override;

    /* make this class non-copyable */
    polynomial_basis(polynomial_basis const&)            = delete;
    polynomial_basis& operator=(polynomial_basis const&) = delete;

    polynomial_basis(polynomial_basis&&) noexcept            = default;
    polynomial_basis& operator=(polynomial_basis&&) noexcept = default;

    uint32_t field_char() const override { return field.n; }

    void clear();

    void clear_redundant()
        requires(not std::same_as<coeff_type, fmpq_class>);

    double memory_usage() const;

    coeff_vect_type const& coefficients(size_t const i) const
    {
        assert(i < m_coefs.size());
        return m_coefs[i];
    }

    void_ptr_type v_coefficients(size_t const i) const override
    {
        assert(i < m_coefs.size());
        return m_coefs[i].data();
    }

    void v_free_generator(size_t const i) override;

#ifdef DEBUG_GAMBA
    void print_generators() const;
#endif

    void export_generators(generators_data& data,
                           bool const only_lead_mons) const;

    void import_generators(generators_data const& data,
                           monomial_order const& mon_order);

    void insert_new_rows_reduce(matrix_f4& matrix,
                                monomial_order const& mon_order);

    void insert_new_rows_echelon(matrix_f4& matrix,
                                 monomial_order const& mon_order,
                                 learn_f4_data* const learn_data,
                                 size_t const round);

    template <class ModularCoefficientType>
    void insert_generators_modular(
        polynomial_basis<ModularCoefficientType>& basis_modp,
        std::vector<std::vector<fmpq_class>>& reconstructed_coefs)
        requires std::same_as<coeff_type, fmpq_class>;

    template <class ModularCoefficientType>
    std::pair<polynomial_basis<ModularCoefficientType>, bool> modular_reduction(
        uint32_t const prime) const;

    void transform_modular_basis()
        requires(not std::same_as<CoefficientType, fmpq_class>);

private:
    void import_generator(size_t const i,
                          generators_data const& data,
                          monomial_order const& mon_order);

    void import_matrix_row(index_vect_type& row,
                           coeff_type* const cfs,
                           std::vector<matrix_monomial_type> const& col_to_mon,
                           monomial_order const& mon_order);

public:
    /* field metadata */
    field_traits<coeff_type> const field;

private:
    /* coefs[i][j] corresponds to the coeff. of the monomial in mons[i][j] */
    std::vector<coeff_vect_type> m_coefs{};
};

template <class CoefficientType>
void polynomial_basis<CoefficientType>::v_free_generator(size_t const i)
{
    assert(i < m_coefs.size());

    /* memory was initialized via a placement new  */
    if constexpr (std::is_same_v<coeff_type, fmpq_class>)
        std::ranges::destroy(m_coefs[i]);

    /* deallocate mons./coefs. memory via the coefficients pointers */
    ::operator delete(m_coefs[i].data());
    m_coefs[i] = {};
}

template <class CoefficientType>
polynomial_basis<CoefficientType>::~polynomial_basis()
{
    for (size_t idx = 0; idx < m_coefs.size(); ++idx)
        v_free_generator(idx);
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::clear()
{
    this->base_polynomial_basis::clear();

    for (size_t idx = 0; idx < m_coefs.size(); ++idx)
        v_free_generator(idx);

    m_coefs.clear();
    m_coefs.shrink_to_fit();
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::clear_redundant()
    requires(not std::same_as<coeff_type, fmpq_class>)
{
    /* copy this before it gets cleared by base class */
    decltype(m_redundant) redundant_copy = m_redundant;

    this->base_polynomial_basis::clear_redundant();

    size_t num_removed = 0;

    for (size_t i = 0; i < redundant_copy.size(); ++i)
    {
        /* clear only non-redundant generators */
        if (not redundant_copy[i])
            continue;

        auto const idx = static_cast<ssize_t>(i - num_removed);

        ::operator delete((std::begin(m_coefs) + idx)->data());
        m_coefs.erase(std::begin(m_coefs) + idx);

        num_removed++;
    }

    /* free memory for deleted redundant generators */
    m_coefs.shrink_to_fit();
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::import_generators(
    generators_data const& data,
    monomial_order const& mon_order)
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
    m_mon_set->reserve(HASHTABLE_INIT_SIZE);

    for (size_t idx = 0; idx < m_num_gens; ++idx)
        import_generator(idx, data, mon_order);

    auto generators       = std::views::zip(m_mons, m_coefs, m_degs);
    auto const gens_begin = std::ranges::begin(generators);

    /* remove empty generators */
    auto const [new_end, end] = std::ranges::remove(
        generators, 0, [](auto const& p) { return get_mons(p).size(); });

    /* deallocate before erasing */
    for (auto const& gen : std::ranges::subrange{new_end, end})
        ::operator delete(get_coef(gen).data());

    ssize_t const new_size = std::distance(gens_begin, new_end);

    /* do the actual removal */
    m_mons.erase(std::cbegin(m_mons) + new_size, std::cend(m_mons));
    m_coefs.erase(std::cbegin(m_coefs) + new_size, std::cend(m_coefs));
    m_degs.erase(std::cbegin(m_degs) + new_size, std::cend(m_degs));

    /* true number of generators after removing empty generators */
    m_num_gens = static_cast<size_t>(new_size);

    /* resize directly data that is constant for all new generators */
    m_redundant.resize(m_num_gens, false);
    m_spair_count.resize(m_num_gens, 0);

    m_num_nondel_gens = m_num_gens;

    m_nondel_gens.resize(m_num_gens);
    /* at the begining no generator is marked as deleted */
    std::iota(std::begin(m_nondel_gens), std::end(m_nondel_gens), 0U);

    generators = std::views::zip(m_mons, m_coefs, m_degs);

    /* sort generators by leading monomial; this ordering will change how the
     * Gebauer-Moeller process creates and updates spairs */
    std::ranges::sort(generators, [&mon_order](auto const& x, auto const& y) {
#if INSERT_ELEMENTS_DECREASING == 1
        return mon_order.cmp(get_mons(x)[0], get_mons(y)[0]) > 0;
#else
        return mon_order.cmp(get_mons(x)[0], get_mons(y)[0]) < 0;
#endif
    });

    /* helper lambda to check if a single polynomials is homogeneous */
    auto const is_homogeneous = [&mon_order](auto const& poly) {
        return std::ranges::all_of(
            get_mons(poly), [&poly, &mon_order](monomial_type const mon) {
                return mon_order.deg(mon) == get_deg(poly);
            });
    };

    /* check if the basis is homogeneous after all the cleaning */
    m_is_homogeneous = std::ranges::all_of(generators, is_homogeneous);

    /* generate divisibility masks after importing generators into basis */
    m_divmap->update_map(*m_mon_set);

    m_lead_sdm.resize(m_num_gens);
    /* compute divmask for each leading term in new basis */
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        monomial_type const lm = m_mons[i][0];
        m_lead_sdm[i]          = m_divmap->compute_divmask(lm);
    }

    m_divmasks_version = m_divmap->version();
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::import_generator(
    size_t const i,
    generators_data const& data,
    monomial_order const& mon_order)
{
    size_t const num_terms = data.lens[i];

    auto const [mon_ptr, cfs_ptr] = allocate_polynomial<coeff_type>(num_terms);

    m_mons[i]  = {mon_ptr, num_terms};
    m_coefs[i] = {cfs_ptr, num_terms};

    // NOLINTNEXTLINE
    size_t const offset = std::accumulate(&data.lens[0], &data.lens[i], 0ULL);

    for (size_t j = 0; j < num_terms; ++j)
    {
        /* convert input monomial into own format & compute hash value */
        auto const mon =
            monomial_init::construct(&data.exps[m_num_vars * (offset + j)],
                                     m_num_vars, params::num_elim_vars);

        /* insert monomial in basis hash table */
        auto const [it, _] = m_mon_set->insert(mon);

        /* copy monomial into generator vector */
        m_mons[i][j] = *it;
    }

    /* import generators coefficients */
    if constexpr (std::is_same_v<coeff_type, fmpq_class>)
    {
        /* allocate_polynomial returns uninitialized memory */
        std::uninitialized_copy(&data.coeffs[offset],
                                &data.coeffs[offset + data.lens[i]],
                                m_coefs[i].data());
    }
    else
    {
        std::copy(&data.coeffs_modp[offset],
                  &data.coeffs_modp[offset + data.lens[i]], m_coefs[i].data());
    }

    auto const poly = std::views::zip(m_mons[i], m_coefs[i]);
    /* sort monomials and coefficients with the given (decreasing) mon. order */
    std::ranges::sort(poly, [&mon_order](auto const x, auto const y) {
        return mon_order.cmp(get_mon(x), get_mon(y)) > 0;
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

    /* remove repeated terms */
    auto const [new_end0, _0] = std::ranges::unique(poly, {}, get_mon);

    /* remove zero terms */
    auto const [new_end, _] =
        std::ranges::remove(poly_begin, new_end0, coeff_type{0}, get_coef);

    /* free mpz integers before shrinking the polynomial's size */
    if constexpr (std::is_same_v<coeff_type, fmpq_class>)
        std::destroy(&get_coef(*new_end), &get_coef(*new_end0));

    /* new size after all simplifications */
    auto const new_size =
        static_cast<size_t>(std::distance(poly_begin, new_end));

    /* do the actual removal by simply updating the size */
    m_mons[i]  = {m_mons[i].data(), new_size};
    m_coefs[i] = {m_coefs[i].data(), new_size};

    /* early return if zero polynomial is found */
    if (new_size == 0)
        return;

    /* normalize coefs by making them monic */
    coeff_type const inv = field.inverse(m_coefs[i][0]);

    if (inv != 1)
    {
        std::ranges::transform(
            m_coefs[i], std::begin(m_coefs[i]),
            [inv, this](coeff_type& c) { return field.multiply(c, inv); });
    }

    assert(m_coefs[i][0] == 1);

    /* transform mod p coefficients to Montgomery space */
    if constexpr (not std::is_same_v<coeff_type, fmpq_class>)
    {
        std::ranges::transform(
            m_coefs[i], std::begin(m_coefs[i]),
            [this](coeff_type const c) { return field.transform(c); });

        assert(m_coefs[i][0] == field.r);
    }

    /* compute degrees after all simplifications */
    if (mon_order.is_degree_order())
    {
        m_degs[i] = mon_order.deg(m_mons[i][0]);
    }
    else
    {
        auto const it = std::ranges::max_element(
            m_mons[i], {}, [&mon_order](monomial_type const mon) {
                return mon_order.deg(mon);
            });

        m_degs[i] = mon_order.deg(*it);
    }

    if (m_degs[i] == 0)
        m_is_trivial = true;
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::insert_new_rows_reduce(
    matrix_f4& matrix,
    monomial_order const& mon_order)
{
    /* timings */
    auto const start_cputime  = std::clock();
    auto const start_walltime = std::chrono::system_clock::now();

    assert(m_num_gens == 0);
    assert(m_mons.size() == 0);
    assert(m_coefs.size() == 0);
    assert(m_redundant.empty());

    m_num_gens = matrix.num_bottom_rows();

    /* all generators are non-redundant & reduced after reduce phase */
    m_redundant.resize(m_num_gens, false);
    m_reduced_gens.resize(m_num_gens);
    std::iota(std::begin(m_reduced_gens), std::end(m_reduced_gens), 0UL);

    /* reserve memory for the new generators */
    m_mons.reserve(m_num_gens);
    m_coefs.reserve(m_num_gens);
    m_degs.reserve(m_num_gens);

    /* rows come presorted from the reduce phase */
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        index_vect_type& row = matrix.m_new_rows[i];
        void_ptr_type vcfs   = matrix.m_new_coefs[i];

        /* basis class must be able to modify coefficients */
        auto* const vcfs_ = const_cast<void*>(vcfs);  // NOLINT
        auto* const cfs   = reinterpret_cast<coeff_type*>(vcfs_);

        import_matrix_row(row, cfs, matrix.m_col_to_mon, mon_order);
    }

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::insert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::insert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::insert_new_rows_echelon(
    matrix_f4& matrix,
    monomial_order const& mon_order,
    learn_f4_data* const learn_data,
    size_t const round)
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
        void_ptr_type vcfs   = matrix.m_new_coefs[i];

        /* basis class must be able to modify coefficients */
        auto* const vcfs_ = const_cast<void*>(vcfs);  // NOLINT
        auto* const cfs   = reinterpret_cast<coeff_type*>(vcfs_);

        import_matrix_row(row, cfs, matrix.m_col_to_mon, mon_order);
    }

    matrix.clear();

    /* once the matrix is clear update deleted generators as generators' memory
     * can be freed safely now */
    update_deleted_gens();

    /* update div. masks after importing new generators into basis */
    m_divmap->update_map(*m_mon_set);

    m_lead_sdm.resize(m_num_gens);

    /* only update sdm of previous generators if new divmap version */
    size_t const start_idx =
        m_divmap->version() != m_divmasks_version ? 0UL : prev_num_gens;

    /* update divmask of leading monomials for all generators since spairs
     * can depend on non-reduced generators */
    for (size_t i = start_idx; i < m_num_gens; ++i)
    {
        if (is_deleted(i))
            continue;

        monomial_type const lm = m_mons[i][0];
        m_lead_sdm[i]          = m_divmap->compute_divmask(lm);
    }

    m_divmasks_version = m_divmap->version();

    /* if learning data is available (i.e. for rational reconstruction), learn
     * data for the new generators from this F4 round */
    learn_new_gens_round(*this, mon_order, num_new_rows, learn_data, round);

    /* timings */
    auto const end_cputime  = std::clock();
    auto const end_walltime = std::chrono::system_clock::now();

    stats::insert_walltime +=
        std::chrono::duration<double>(end_walltime - start_walltime).count();
    stats::insert_cputime +=
        static_cast<double>(end_cputime - start_cputime) / CLOCKS_PER_SEC;
}

template <class CoefficientType>
template <class ModularCoefficientType>
void polynomial_basis<CoefficientType>::insert_generators_modular(
    polynomial_basis<ModularCoefficientType>& basis_modp,
    std::vector<std::vector<fmpq_class>>& reconstructed_coefs)
    requires std::same_as<coeff_type, fmpq_class>
{
    assert(m_num_gens == 0);
    assert(m_mons.size() == 0);
    assert(m_coefs.size() == 0);
    assert(m_redundant.empty());

    m_num_gens = basis_modp.m_reduced_gens.size();

    /* all generators are non-redundant & reduced after reduce phase */
    m_redundant.resize(m_num_gens, false);
    m_reduced_gens.resize(m_num_gens);
    std::iota(std::begin(m_reduced_gens), std::end(m_reduced_gens), 0UL);

    /* resize memory for reconstructed generators */
    m_mons.resize(m_num_gens);
    m_coefs.resize(m_num_gens);
    m_degs.resize(m_num_gens);

    for (size_t i = 0; i < basis_modp.m_reduced_gens.size(); ++i)
    {
        size_t const idx = basis_modp.m_reduced_gens[i];
        size_t const len = basis_modp.length(idx);

        assert(reconstructed_coefs[idx].size() != 0);
        assert(reconstructed_coefs[idx].size() == len);

        auto const [mon_ptr, cfs_ptr] = allocate_polynomial<coeff_type>(len);

        monomial_vect_type const& mons = basis_modp.m_mons[idx];
        /* copy basis_modp monomials to final rational basis */
        std::copy(std::begin(mons), std::end(mons), mon_ptr);

        std::vector<fmpq_class>& rec_coefs = reconstructed_coefs[idx];
        /* allocate_polynomial returns uninitialized memory; move reconstructed
         * coefficients to final basis over the rationals */
        std::uninitialized_copy(std::make_move_iterator(std::begin(rec_coefs)),
                                std::make_move_iterator(std::end(rec_coefs)),
                                cfs_ptr);

        m_mons[i]  = {mon_ptr, len};
        m_coefs[i] = {cfs_ptr, len};
        m_degs[i]  = basis_modp.m_degs[idx];
    }
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::import_matrix_row(
    polynomial_basis<CoefficientType>::index_vect_type& row,
    polynomial_basis<CoefficientType>::coeff_type* const cfs,
    std::vector<matrix_monomial_type> const& col_to_mon,
    monomial_order const& mon_order)
{
    size_t const num_terms = row.size();

    auto const* const ind_ptr = reinterpret_cast<index_ptr_type>(row.data());
    auto* const mon_ptr       = reinterpret_cast<monomial_type*>(row.data());

    monomial_vect_type mons = {mon_ptr, num_terms};

    for (size_t i = 0; i < num_terms; ++i)
    {
        assert(ind_ptr[i] < col_to_mon.size());
        /* undo the columns <-> monomials transformation */
        matrix_monomial_type mat_mon = col_to_mon[ind_ptr[i]];

        /* copy matrix monomial into a basis monomial */
        auto const mon = monomial_init::copy(mat_mon);

        /* insert basis monomial into basis hash table */
        auto const [it, _] = m_mon_set->insert(mon);

        /* copy monomial to generator vector */
        mons[i] = *it;
    }

    m_mons.emplace_back(mons);
    m_coefs.emplace_back(cfs, num_terms);

    /* compute degree of new generator */
    if (mon_order.is_degree_order())
    {
        size_t const deg = mon_order.deg(m_mons.back()[0]);

        m_degs.emplace_back(deg);
    }
    else
    {
        auto const it = std::ranges::max_element(
            m_mons.back(), {}, [&mon_order](monomial_type const mon) {
                return mon_order.deg(mon);
            });

        size_t const deg = mon_order.deg(*it);

        m_degs.emplace_back(deg);
    }

    if (m_degs.back() == 0)
        m_is_trivial = true;
}

template <>
template <class ModularCoefficientType>
std::pair<polynomial_basis<ModularCoefficientType>, bool>
polynomial_basis<fmpq_class>::modular_reduction(uint32_t const prime) const
{
    using mod_coeff_type      = ModularCoefficientType;
    using mod_basis_type      = polynomial_basis<mod_coeff_type>;
    using mod_coeff_vect_type = mod_basis_type::coeff_vect_type;
    /* make sure prime number fits in the requested coefficient type */
    assert(prime < std::numeric_limits<mod_coeff_type>::max());
    /* absolute maximum supported field characteristic is 2^31 - 1 */
    assert(prime <= std::numeric_limits<int32_t>::max());

    /* first update the lead_sdm if 'divmap' has been updated by other bases */
    if (m_divmasks_version != m_divmap->version())
    {
        for (size_t i = 0; i < m_num_gens; ++i)
        {
            monomial_type const lm = m_mons[i][0];
            m_lead_sdm[i]          = m_divmap->compute_divmask(lm);
        }

        m_divmasks_version = m_divmap->version();
    }

    /* partially construct new basis modulo the given prime */
    polynomial_basis<mod_coeff_type> basis_modp{prime};
    field_traits<mod_coeff_type> const& field_modp = basis_modp.field;

    /* set base_polynomial_basis member variables */
    basis_modp.m_num_vars = m_num_vars;
    basis_modp.m_num_gens = m_num_gens;
    basis_modp.m_mon_set  = m_mon_set;
    // basis_modp.m_mons     = m_mons; // don't share monomials
    basis_modp.m_degs = m_degs;

    basis_modp.m_divmap           = m_divmap;
    basis_modp.m_lead_sdm         = m_lead_sdm;
    basis_modp.m_divmasks_version = m_divmasks_version;

    basis_modp.m_redundant       = m_redundant;
    basis_modp.m_reduced_gens    = m_reduced_gens;
    basis_modp.m_num_nondel_gens = m_num_gens;
    basis_modp.m_spair_count     = m_spair_count;
    basis_modp.m_nondel_gens     = m_nondel_gens;

    /* the definition of bad prime implies these flags don't change */
    basis_modp.m_is_homogeneous = m_is_homogeneous;
    basis_modp.m_is_trivial     = m_is_trivial;

    bool is_bad_prime = false;

    basis_modp.m_coefs.reserve(m_num_gens);
    /* reduce coefs. modulo the given prime & move them to Montgomery space */
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        coeff_vect_type const& cfs     = m_coefs[i];
        monomial_vect_type const& mons = m_mons[i];

        size_t const num_terms = cfs.size();

        auto const [mon_ptr, cfs_ptr] =
            allocate_polynomial<mod_coeff_type>(num_terms);

        monomial_vect_type& mons_modp =
            basis_modp.m_mons.emplace_back(mon_ptr, num_terms);
        mod_coeff_vect_type& cfs_modp =
            basis_modp.m_coefs.emplace_back(cfs_ptr, num_terms);

        /* copy monomials to modular basis */
        std::copy(std::begin(mons), std::end(mons), std::begin(mons_modp));

        /* do the actual reduction modulo the given prime */
        std::ranges::transform(
            cfs, std::begin(cfs_modp), [&field_modp](fmpq_class const& c) {
                uint32_t const a = field_modp.modular_reduce(c.get_num());
                uint32_t const b = field_modp.modular_reduce(c.get_den());
                return field_modp.multiply(a, field_modp.inverse(b));
            });

        /* modular reduced polynomials are already monic */
        assert(cfs_modp[0] == 1);

        /* a prime is bad if it divides any coefficient of the bais */
        is_bad_prime = std::ranges::any_of(
            cfs_modp, [](mod_coeff_type const c) { return c == 0; });

        /* early break, return a partially constructed basis */
        if (is_bad_prime)
            break;

        /* transform mod p coefficients to Montgomery space */
        std::ranges::transform(cfs_modp, std::begin(cfs_modp),
                               [&field_modp](mod_coeff_type const c) {
                                   return field_modp.transform(c);
                               });

        assert(cfs_modp[0] == basis_modp.field.r);
    }

    return std::make_pair(std::move(basis_modp), is_bad_prime);
}

#ifdef DEBUG_GAMBA
template <class CoefficientType>
void polynomial_basis<CoefficientType>::print_generators() const
{
    for (size_t idx = 0; idx < m_num_gens; ++idx)
    {
        std::string gen_str;
        gen_str += fmt::format("#{}: deg = {}, ", idx, m_degs[idx]);

        for (size_t j = 0; j < m_mons[idx].size(); ++j)
        {
            std::string const& mon_str =
                monomial2string<monomial_type>(m_mons[idx][j]);

            gen_str += (j != 0 ? " + " : "");

            if constexpr (std::is_same_v<coeff_type, fmpq_class>)
                gen_str += m_coefs[idx][j].get_str();
            else
                gen_str += fmt::format("{:d}", m_coefs[idx][j]);

            gen_str += (mon_str.empty() ? "" : "*") + mon_str;
        }

        log::print(log::DEBG, "{}\n", gen_str);
    }
}
#endif

template <class CoefficientType>
void polynomial_basis<CoefficientType>::export_generators(
    generators_data& data,
    bool const only_lead_mons) const
{
    using var_type     = generators_data::var_type;
    using coeff_p_type = generators_data::coeff_p_type;

    size_t const num_red_gens = m_reduced_gens.size();

    data.num_vars   = static_cast<var_type>(m_num_vars);
    data.field_char = static_cast<var_type>(m_field_char);
    data.num_gens   = static_cast<var_type>(num_red_gens);

    size_t const num_total_mons =
        only_lead_mons ? m_mons.size()
                       : std::accumulate(std::cbegin(m_mons), std::cend(m_mons),
                                         0ULL, [](size_t acc, auto const& v) {
                                             return acc + v.size();
                                         });

    if constexpr (std::is_same_v<coeff_type, fmpq_class>)
    {
        data.coeffs.resize(num_total_mons);
    }
    else
    {
        data.coeffs_modp.resize(num_total_mons);
    }

    data.exps.resize(num_total_mons * m_num_vars);
    data.lens.resize(data.num_gens);

    /* write generators in increasing order of leading monomial */
    for (size_t i = 0, offset = 0; i < num_red_gens; ++i)
    {
        size_t const idx = m_reduced_gens[i];

        data.lens[i] = only_lead_mons ? 1UL : m_mons[idx].size();

        for (size_t j = 0; j < data.lens[i]; ++j)
        {
            monomial2exponent(m_mons[idx][j],
                              &data.exps[(offset + j) * m_num_vars]);

            if constexpr (std::is_same_v<coeff_type, fmpq_class>)
            {
                data.coeffs[offset + j] = m_coefs[idx][j];
            }
            else
            {
                /* transform mod p coefficients from Montgomery space */
                coeff_type const c = field.reduce_normalize(m_coefs[idx][j]);

                data.coeffs_modp[offset + j] = static_cast<coeff_p_type>(c);
            }
        }

        offset += data.lens[i];
    }
}

template <class CoefficientType>
void polynomial_basis<CoefficientType>::transform_modular_basis()
    requires(not std::same_as<CoefficientType, fmpq_class>)
{
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        /* transform mod p coefficients out of Montgomery space */
        for (size_t j = 0; j < m_mons[i].size(); ++j)
            m_coefs[i][j] = field.reduce_normalize(m_coefs[i][j]);
    }
}

template <class CoefficientType>
double polynomial_basis<CoefficientType>::memory_usage() const
{
    double mem_size = base_polynomial_basis::memory_usage();

    mem_size += memory_size(m_coefs);
    for (size_t i = 0; i < m_num_gens; ++i)
    {
        mem_size += memory_size(m_coefs[i]);
    }

    return mem_size;
}

}  // namespace gamba
