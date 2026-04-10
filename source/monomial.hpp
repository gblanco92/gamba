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
#include <cassert>
#include <random>

#include "divmask.hpp"
#include "params.hpp"
#include "utils.hpp"

namespace gamba
{

/* monomial context classes */
struct basis_hashtable
{};

struct spair_hashtable
{};

struct matrix_hashtable
{};

template <class T>
struct monomial_data
{
    using degree_type = uint32_t;
    using hash_type   = uint32_t;
    using index_type  = uint32_t;
    using count_type  = uint32_t;
};

template <>
struct monomial_data<basis_hashtable> : monomial_data<void>
{
    hash_type hash;
}; /* end struct monomial_data for basis_hashtable */

template <>
struct monomial_data<spair_hashtable> : monomial_data<void>
{
    hash_type hash;
}; /* end struct monomial_data for spair_hashtable */

template <>
struct monomial_data<matrix_hashtable> : monomial_data<void>
{
    hash_type hash;
    index_type idx;
}; /* end struct monomial_data for matrix_hashtable */

struct monomial_base
{
    using exponent_type     = uint16_t;
    using exponent_ptr_type = exponent_type const*;

    using monomial_data_type = monomial_data<void>;
    using hash_type          = monomial_data_type::hash_type;
    using weight_type        = hash_type;

    size_t size() const { return exp_size; }  // NOLINT

    static void initialize_weights()
    {
        std::mt19937_64 eng{params::seed};

        std::uniform_int_distribution<weight_type> dist(
            1, std::numeric_limits<weight_type>::max());

        m_weights.resize(monomial_base::exp_size);
        std::generate(std::begin(m_weights), std::end(m_weights),
                      [&]() { return dist(eng); });
    }

    /* length of exponent vector including degree(s) */
    static inline size_t exp_size{0};
    /* number of elimination variables + degree of the block */
    static inline size_t block_size{0};

protected:
    /* weights to compute the hashes */
    static inline std::vector<weight_type> m_weights;
}; /* end class monomial_base */

/* forward declaration */
template <class>
class monomial;

template <class MonomialContext>
class monomial_init : public monomial_base
{
public:
    using monomial_context = MonomialContext;
    using monomial_type    = monomial<monomial_context>;

private:
    /* constructor only accessible via static methods */
    monomial_init() { assert(m_exp != nullptr); }

public:
    exponent_ptr_type cbegin() const { return m_exp; }

    exponent_ptr_type cend() const { return m_exp + size(); }

    /* the monomial_init template class has always its hash value stored in the
     * member variable hash */
    hash_type hash() const { return m_hash; }

private:
    FORCE_INLINE void compute_hash()
    {
        m_hash = 0;
        for (size_t i = 0; i < m_weights.size(); ++i)
            m_hash += m_weights[i] * m_exp[i];
    }

public:
    /* construct a valid monomial from external memory */
    static monomial_init construct(int32_t const* const mon_exp,
                                   size_t const num_vars,
                                   size_t const num_elim_vars)
    {
        size_t const off = num_elim_vars > 0 ? 2 : 1;
        size_t const ebk = block_size;

        monomial_init mon;
        m_exp[0] = m_exp[ebk] = 0;

        for (size_t i = 0; i < num_elim_vars; ++i)
        {
            m_exp[i + 1] = static_cast<exponent_type>(mon_exp[i]);
            /* degree */
            m_exp[0] += m_exp[i + 1];
        }
        for (size_t i = num_elim_vars; i < num_vars; ++i)
        {
            m_exp[i + off] = static_cast<exponent_type>(mon_exp[i]);
            /* degree */
            m_exp[ebk] += m_exp[i + off];
        }

        mon.compute_hash();

        return mon;
    }

    /* creates the least common multiple (lcm) of the lhs and rhs monomials */
    template <class OtherMonomialContext>
    static monomial_init lcm(monomial<OtherMonomialContext> const lhs,
                             monomial<OtherMonomialContext> const rhs)
        requires std::same_as<monomial_context, spair_hashtable>
    {
        size_t const esz = exp_size;
        size_t const ebk = block_size;

        exponent_ptr_type const elhs = lhs.cbegin();
        exponent_ptr_type const erhs = rhs.cbegin();

        monomial_init mon;
        /* ignore degree(s) */
        for (size_t i = 1; i < esz; ++i)
            m_exp[i] = elhs[i] < erhs[i] ? erhs[i] : elhs[i];

        /* reset degrees */
        m_exp[0] = m_exp[ebk] = 0;

        /* if not using an elimination order the first loop is not executed and
         * the degree is computed in the second loop */
        for (size_t i = 1; i < ebk; ++i)
            m_exp[0] += m_exp[i];
        for (size_t i = ebk + 1; i < esz; ++i)
            m_exp[ebk] += m_exp[i];

        mon.compute_hash();

        return mon;
    }

    /* creates the quotient monomial mon1/mon2 assuming mon2 divides mon1 */
    template <class OtherMonomialContext>
    static monomial_init quotient(monomial<OtherMonomialContext> const mon1,
                                  monomial_type const mon2)
    {
        size_t const esz = exp_size;

        exponent_ptr_type const exp1 = mon1.cbegin();
        exponent_ptr_type const exp2 = mon2.cbegin();

        monomial_init mon;
        mon.m_hash = mon1.hash() - mon2.hash();

        for (size_t i = 0; i < esz; ++i)
        {
            assert(exp1[i] >= exp2[i]);
            m_exp[i] = exp1[i] - exp2[i];
        }

        return mon;
    }

    /* create the quotient monomial mon1 * mon2, may overflow exponent type */
    template <class OtherMonomialContext>
    static monomial_init product(monomial_init<basis_hashtable> const mon1,
                                 monomial<OtherMonomialContext> const mon2)
        requires std::same_as<monomial_context, matrix_hashtable>
    {
        size_t const esz = exp_size;

        exponent_ptr_type const exp1 = mon1.cbegin();
        exponent_ptr_type const exp2 = mon2.cbegin();

        monomial_init mon;
        mon.m_hash = mon1.hash() + mon2.hash();

        for (size_t i = 0; i < esz; ++i)
            m_exp[i] = exp1[i] + exp2[i];

        return mon;
    }

    template <class OtherMonomialContext>
    static monomial_init copy(monomial<OtherMonomialContext> const other)
        requires(not std::same_as<monomial_context, OtherMonomialContext>)
    {
        size_t const esz = exp_size;

        exponent_ptr_type const exp = other.cbegin();

        monomial_init mon;
        /* all monomial contexts share the same weights */
        mon.m_hash = other.hash();

        for (size_t i = 0; i < esz; ++i)
            m_exp[i] = exp[i];

        return mon;
    }

private:
    template <class>
    friend class monomial;

    template <class>
    friend class monomial_flat_set;

    /* the exponent table is local to each monomial context */
    static inline exponent_type* m_exp{nullptr};  // TODO thread_local

    hash_type m_hash{};
}; /* end class monomial_init */

template <class MonomialContext>
class monomial : public monomial_base
{
public:
    using monomial_context = MonomialContext;
    using monomial_data_t  = monomial_data<monomial_context>;
    using monomial_init_t  = monomial_init<monomial_context>;

    using index_type  = typename monomial_data_t::index_type;
    using hash_type   = typename monomial_data_t::hash_type;
    using degree_type = typename monomial_data_t::degree_type;
    using count_type  = typename monomial_data_t::count_type;

    monomial() noexcept : m_idx{} {}

    /* explicitly constructor for hash table insertion; store the hash value
     * coming from the monomial_init so that the hash value is stored in-place
     * in the hash table */
    explicit monomial(monomial_init_t const mon) :
            m_idx{static_cast<index_type>(load)}
    {
        data_v[m_idx].hash = mon.hash();

        if constexpr (std::is_same_v<monomial_context, matrix_hashtable>)
        {
            data_v[m_idx].idx = 0;
            // data_v[m_idx].cnt = 0;
        }

        /* the hash table always reserves an extra spot so exp always points
         * to a valid memory location */
        monomial_init_t::m_exp += size();  // TODO thread_local
        monomial::load++;
    }

    static constexpr monomial dummy()
    {
        monomial dummy;
        dummy.m_idx = std::numeric_limits<index_type>::max();

        return dummy;
    }

    bool operator==(monomial const other) const { return m_idx == other.m_idx; }

    exponent_ptr_type cbegin() const { return exps_v + size() * m_idx; }

    exponent_ptr_type cend() const { return exps_v + size() * (m_idx + 1); }

    hash_type hash() const { return data_v[m_idx].hash; }

    monomial_data_t& data() const { return data_v[m_idx]; }

    /* returns true if lhs divides rhs (lhs | rhs), otherwise returns false */
    template <class OtherMonomialContext>
    static bool check_monomial_division(
        monomial<OtherMonomialContext> const lhs,
        monomial const rhs)
    {
        exponent_ptr_type const elhs = lhs.cbegin();
        exponent_ptr_type const erhs = rhs.cbegin();

        exponent_type flag{0};
        for (size_t i = 0; i < exp_size; ++i)
            flag |= (elhs[i] > erhs[i]);

        return flag == 0;
    }

    /* returns true if lhs divides rhs (lhs | rhs), otherwise returns false */
    template <class OtherMonomialContext>
    static bool check_monomial_division_divmask(
        monomial<OtherMonomialContext> lhs_mon,
        divmask_type const lhs_sdm,
        monomial const rhs_mon,
        divmask_type const rhs_sdm)
    {
        /* use divmasks for fast divisibility check */
        if (not(lhs_sdm <= rhs_sdm))
        {
            assert(check_monomial_division(lhs_mon, rhs_mon) == 0);
            return false;
        }

        exponent_ptr_type const elhs = lhs_mon.cbegin();
        exponent_ptr_type const erhs = rhs_mon.cbegin();

        exponent_type flag{0};
        for (size_t i = 0; i < exp_size; ++i)
            flag |= (elhs[i] > erhs[i]);

        return flag == 0;
    }

    /* returns true if gcd(lhs, rhs) = 1 */
    static bool coprime_monomials(monomial const lhs, monomial const rhs)
    {
        exponent_ptr_type const elhs = lhs.cbegin();
        exponent_ptr_type const erhs = rhs.cbegin();

        exponent_type flag{0};
        /* take into account the presence of degree(s) */
        for (size_t i = 1; i < block_size; ++i)
            flag |= (elhs[i] != 0 and erhs[i] != 0);
        for (size_t i = block_size + 1; i < exp_size; ++i)
            flag |= (elhs[i] != 0 and erhs[i] != 0);

        return flag == 0;
    }

    static exponent_ptr_type exps_vect() { return exps_v; }

private:
    template <class>
    friend class monomial_flat_set;

    static inline monomial_data_t* data_v{nullptr};
    static inline exponent_type* exps_v{nullptr};

    static inline size_t load;
    // static inline size_t idx; // TODO thread_local

    index_type m_idx;
}; /* end class monomial */

struct monomial_equal_to
{
    using exponent_type     = monomial_base::exponent_type;
    using exponent_ptr_type = monomial_base::exponent_ptr_type;

    using is_transparent = void;

    template <class MonomialContext>
    bool operator()(monomial_init<MonomialContext> const mon1,
                    monomial<MonomialContext> const mon2) const
    {
        if (mon1.hash() != mon2.hash())
            return false;

        exponent_ptr_type const exp1 = mon1.cbegin();
        exponent_ptr_type const exp2 = mon2.cbegin();

        exponent_type flag{0};
        for (size_t i = 0; i < mon1.size(); ++i)
            flag |= (exp1[i] ^ exp2[i]);

        return flag == 0;
    }
}; /* end class monomial_equal_to */

}  // namespace gamba
