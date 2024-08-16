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

#include <bit>
#include <memory>

#include "allocator.hpp"
#include "config.hpp"
#include "monomial.hpp"
#include "stats.hpp"
#include "utils.hpp"

namespace gamba
{

template <class T>
class monomial_flat_set
{};

template <class MonomialContext, class MonomialOrder>
class monomial_flat_set<monomial<MonomialContext, MonomialOrder>>
{
    static constexpr double const m_max_load_factor = HASH_LOAD_FACTOR;
    static constexpr size_t const m_inv_max_load_factor =
        static_cast<size_t>(1.0 / m_max_load_factor);

    static_assert(std::has_single_bit(m_inv_max_load_factor),
                  "Inverse of max. load factor must be a power of 2.");

public:
    using monomial_context = MonomialContext;
    using monomial_order   = MonomialOrder;
    using monomial_type    = monomial<monomial_context, monomial_order>;
    using monomial_data    = typename monomial_type::monomial_data_t;
    using monomial_init    = typename monomial_type::monomial_init_t;
    using index_type       = typename monomial_type::index_type;
    using exponent_type    = typename monomial_type::exponent_type;

private:
    class iterator
    {
    public:
        using difference_type   = std::ptrdiff_t;
        using value_type        = monomial_type;
        using pointer           = monomial_type const*;
        using reference         = monomial_type const&;
        using iterator_category = std::bidirectional_iterator_tag;

        explicit iterator(size_t const i) : m_idx{static_cast<index_type>(i)} {}

        explicit iterator(monomial_type const mon) :
                m_idx{*reinterpret_cast<index_type const*>(&mon)}
        {}

        iterator& operator++()
        {
            ++m_idx;
            return *this;
        }

        iterator operator++(int)
        {
            iterator retval = *this;
            ++(*this);
            return retval;
        }

        iterator operator+(size_t const delta) const
        {
            iterator retval = *this;
            retval.m_idx += static_cast<index_type>(delta);
            return retval;
        }

        iterator& operator--()
        {
            --m_idx;
            return *this;
        }

        iterator operator--(int)
        {
            iterator retval = *this;
            --(*this);
            return retval;
        }

        iterator operator-(size_t const delta) const
        {
            iterator retval = *this;
            retval.m_idx -= static_cast<index_type>(delta);
            return retval;
        }

        bool operator==(iterator const other) const
        {
            return m_idx == other.m_idx;
        }

        bool operator!=(iterator const other) const
        {
            return !(*this == other);
        }

        pointer operator->() const { return reinterpret_cast<pointer>(&m_idx); }

        reference operator*() const { return *operator->(); }

    private:
        index_type m_idx;
    }; /* end class iterator */

public:
    using key_type        = monomial_type;
    using value_type      = key_type;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;
    // using hasher          = monomial_hash;
    using key_equal = monomial_equal_to;

    static constexpr size_t const huge_page_size =
        aligned_allocator<int>::huge_page_size;
    using allocator_type = aligned_allocator<monomial_type, huge_page_size>;

    using reference       = value_type&;
    using const_reference = value_type const&;
    using pointer = typename std::allocator_traits<allocator_type>::pointer;
    using const_pointer =
        typename std::allocator_traits<allocator_type>::const_pointer;
    using iterator       = monomial_flat_set::iterator;
    using const_iterator = monomial_flat_set::iterator;

    using allocator_data_type = std::allocator_traits<
        allocator_type>::template rebind_alloc<monomial_data>;
    using allocator_exps_type = std::allocator_traits<
        allocator_type>::template rebind_alloc<exponent_type>;

    monomial_flat_set()
    {
        /* make sure it is always safe to construct one object */
        reserve(1);

        monomial_type::load = 0;
    }

    /* make this class non-copyable */
    monomial_flat_set(monomial_flat_set const&)            = delete;
    monomial_flat_set& operator=(monomial_flat_set const&) = delete;

    void clear()
    {
        m_size = 0;

        std::fill(std::assume_aligned<CACHE_LINE_SIZE>(m_table),
                  std::assume_aligned<CACHE_LINE_SIZE>(m_table + m_capacity),
                  monomial_type::dummy());

        monomial_type::load  = 0;
        monomial_init::m_exp = &monomial_type::exps_v[0];  // TODO thread_local
    }

    ~monomial_flat_set()
    {
        m_size = 0;

        m_alloc.deallocate(m_table, m_capacity);
        m_table = nullptr;

        m_alloc_data.deallocate(monomial_type::data_v, m_max_size);
        monomial_type::data_v = nullptr;

        m_alloc_exps.deallocate(monomial_type::exps_v,
                                monomial_base::exp_size * m_max_size);
        monomial_type::exps_v = nullptr;

        monomial_type::load = 0;
    }

    size_t size() const
    {
        assert(m_size == monomial_type::load);
        return m_size;
    }

    size_t max_size() const { return m_max_size; }

    void rehash(size_t const capacity)
    {
        /* since we are going to reserve the whole page anyway, make hash table
         * capacity *in bytes* a multiple of the huge page size */
        size_t const huge_capacity =
            round_up(sizeof(key_type) * capacity, huge_page_size)
            / sizeof(key_type);

        /* space for the exponents according to the max_load_factor */
        if (huge_capacity > m_capacity)
        {
            rehash_impl(huge_capacity, huge_capacity / m_inv_max_load_factor);

            if constexpr (std::is_same_v<monomial_context, basis_hashtable>)
                stats.max_size_bht = std::max(stats.max_size_bht, m_max_size);

            if constexpr (std::is_same_v<monomial_context, spair_hashtable>)
                stats.max_size_sht = std::max(stats.max_size_sht, m_max_size);

            if constexpr (std::is_same_v<monomial_context, matrix_hashtable>)
                stats.max_size_mht = std::max(stats.max_size_mht, m_max_size);
        }
    }

    void reserve(size_t const count)
    {
        /* reserve enough space for hashing at least 'count' num. monomials */
        rehash(count * m_inv_max_load_factor);
    }

    void prefetch_insert(size_t const hash) const
    {
        size_t const k = hash & (m_capacity - 1);
        _mm_prefetch(m_table + k, _MM_HINT_T0);
    }

    std::pair<iterator, bool> insert(monomial_init const mon)
    {
        /* capacity is always a power of two */
        size_t const mod = m_capacity - 1;

        size_t k = mon.hash(), i = 0;
        while (true)
        {
            /* linear proving */
            k = (k + i++) & mod;

            monomial_type const mon_k = m_table[k];

            /* empty position */
            if (mon_k == monomial_type::dummy())
                break;

            /* monomial already in hash table */
            if (LIKELY(monomial_equal_to{}(mon, mon_k)))
                return {iterator{mon_k}, false};
        }

        /* insert new element */
        monomial_type const new_mon = (m_table[k] = monomial_type{mon});

        m_size++;

        /* rehash here since we have no place for creating temp. exponents */
        if (UNLIKELY(m_size == m_max_size))
            reserve(2 * m_size);

        /* return new_mon that is not invalidated after rehashing */
        return {iterator{new_mon}, true};
    }

    iterator begin() const { return iterator{0}; }

    iterator end() const { return iterator{monomial_type::load}; }

    iterator cbegin() const { return iterator{0}; }

    iterator cend() const { return iterator{monomial_type::load}; }

    double load_factor() const
    {
        return static_cast<double>(m_size) / static_cast<double>(m_capacity);
    }

    static constexpr double max_load_factor() { return m_max_load_factor; }

    double memory_usage() const
    {
        double mem_size = 0;

        mem_size += static_cast<double>(m_capacity * sizeof(index_type))
                  / 1024.0 / 1024.0;
        mem_size += static_cast<double>(m_max_size * sizeof(monomial_data))
                  / 1024.0 / 1024.0;
        mem_size += static_cast<double>(m_max_size * monomial_base::exp_size
                                        * sizeof(exponent_type))
                  / 1024.0 / 1024.0;

        return mem_size;
    }

private:
    void rehash_impl(size_t const capacity_hash, size_t const capacity_exp)
    {
        /* capacity_hash must be a power of two for fast modulo */
        assert(std::has_single_bit(capacity_hash));
        /* capacity_exp must be a power of two for packing routines */
        assert(std::has_single_bit(capacity_exp));

        m_alloc.deallocate(m_table, m_capacity);
        m_table    = m_alloc.allocate(capacity_hash);
        m_capacity = capacity_hash;

        std::fill(std::assume_aligned<CACHE_LINE_SIZE>(m_table),
                  std::assume_aligned<CACHE_LINE_SIZE>(m_table + m_capacity),
                  monomial_type::dummy());

        size_t const mod = m_capacity - 1;

        for (auto it = cbegin(); it != cend(); ++it)
        {
            size_t k = it->hash(), j = 0;
            while (true)
            {
                /* linear proving */
                k = (k + j++) & mod;

                if (m_table[k] != monomial_type::dummy())
                    continue;

                m_table[k] = *it;

                break;
            }
        }

        /* allocate and reallocate data */
        monomial_data* const new_data_v = m_alloc_data.allocate(capacity_exp);

        std::copy(std::assume_aligned<CACHE_LINE_SIZE>(monomial_type::data_v),
                  std::assume_aligned<CACHE_LINE_SIZE>(monomial_type::data_v
                                                       + monomial_type::load),
                  new_data_v);

        m_alloc_data.deallocate(monomial_type::data_v, m_max_size);
        monomial_type::data_v = new_data_v;

        /* allocate and reallocate exps */
        exponent_type* const new_exps_v =
            m_alloc_exps.allocate(capacity_exp * monomial_base::exp_size);

        std::copy(std::assume_aligned<CACHE_LINE_SIZE>(monomial_type::exps_v),
                  std::assume_aligned<CACHE_LINE_SIZE>(
                      monomial_type::exps_v
                      + monomial_type::load * monomial_base::exp_size),
                  new_exps_v);

        m_alloc_exps.deallocate(monomial_type::exps_v,
                                m_max_size * monomial_base::exp_size);
        monomial_type::exps_v = new_exps_v;

        /* update variables */
        m_max_size = capacity_exp;

        // TODO: thread_local
        monomial_init::m_exp = monomial_type::exps_v
                             + (monomial_type::load * monomial_base::exp_size);
    }

    /* allocators */
    allocator_type m_alloc{};
    allocator_data_type m_alloc_data{};
    allocator_exps_type m_alloc_exps{};

    /* the hash space */
    key_type* m_table{nullptr};

    /* the capacity of the hash table; must be a power of two */
    size_t m_capacity{0};

    /* the maximum number of elements before we need to rehash */
    size_t m_max_size{0};

    /* number of elements currently present in the hash table */
    size_t m_size{0};
}; /* end monomial_flat_set class */

}  // namespace gamba
