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

#include <cstring>
#include <ranges>
#include <string>

#include <gmpxx.h>

#include "monomial.hpp"

namespace gamba
{

struct generators_data
{
    using var_type     = uint32_t;
    using coeff_type   = mpz_class;
    using coeff_p_type = uint32_t;
    using exp_type     = int32_t;

    /* type used by the engine to represent exponent */
    using exponent_type = typename monomial_base::exponent_type;

    void read(std::istream& infile);

private:
    void read_num_vars(std::istream& infile);

    void read_characteristic(std::istream& infile);

    void read_num_generators(std::istream& infile);

    void read_variable_names(std::istream& infile);

    void read_generators(std::istream& infile);

    void read_generator_line(std::string line, size_t const line_num);

    void read_exponent(std::string term);

public:
    void write(std::ostream& outfile) const;

private:
    template <class CoefficientType>
    void write_generators(std::ostream& outfile,
                          std::vector<CoefficientType> const& coeff) const;
    template <bool Star>
    bool write_monomial(std::ostream& outfile, size_t offset) const;

public:
    var_type num_vars{};
    var_type field_char{};
    var_type num_gens{};

    std::vector<std::string> var_names{};
    std::vector<std::string> sorted_var_names{};
    std::vector<size_t> argsort_vnames{};

    std::vector<coeff_type> coeffs{};
    std::vector<coeff_p_type> coeffs_modp{};
    std::vector<exp_type> exps{};
    std::vector<size_t> lens{};
};

template <class MonomialType>
std::string monomial2string(MonomialType const& mon)
{
    using monomial_type  = MonomialType;
    using monomial_order = typename monomial_type::monomial_order;

    constexpr size_t const offset =
        std::is_same_v<monomial_order, order_blockelim> ? 2 : 1;

    static auto const names =
        std::views::iota(0ULL, monomial_type::size() - offset)
        | std::views::transform([](auto i) { return "x" + std::to_string(i); });

    auto const* const exp = mon.cbegin();

    std::string str;
    bool first{true};

    for (size_t i = 1; auto name : names)
    {
        auto const e = exp[i++];
        /* skip degree of second block (if present) */
        if (i == monomial_order::block_size)
            ++i;

        if (e == 0)
            continue;

        str += (first ? "" : "*") + name + "^" + std::to_string(e);

        first = false;
    }

    return str.empty() ? "1" : str;
}

template <class MonomialType>
void monomial2exponent(MonomialType const& mon,
                       generators_data::exp_type* const out)
{
    using monomial_type  = MonomialType;
    using monomial_order = typename monomial_type::monomial_order;

    auto const* const exp = mon.cbegin();

    for (size_t i = 1, j = 0; i < mon.size(); ++i)
    {
        /* skip degree of second block (if present) */
        if (i == monomial_order::block_size)
            continue;

        out[j++] = static_cast<generators_data::exp_type>(exp[i]);
    }
}

}  // namespace gamba
