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

#include "io.hpp"

#include <flint/ulong_extras.h>

#include "thirdparty/flintxx.hpp"

namespace gamba
{

namespace
{

void remove_spaces(std::string& line)
{
    auto is_space = [](unsigned char c) { return std::isspace(c); };

    line.erase(std::remove_if(std::begin(line), std::end(line), is_space),
               std::cend(line));
}

void remove_trailing_comma(std::string& line)
{
    auto const pos = line.find_last_of(',');

    if (pos != std::string::npos and pos == line.size() - 1)
        line.resize(line.size() - 1);
}

std::string read_term(std::string& line, std::string::iterator& it)
{
    static std::string const str{"+-"};
    auto next = std::find_first_of(it + 1, std::end(line), std::cbegin(str),
                                   std::cend(str));

    std::string term{it, next};
    it = next;
    return term;
}

fmpq_class term_to_mpq(std::string line, size_t* pos)
{
    /* skip and save sign in case the monomial is -x or +y */
    char sign{'\0'};
    if (line[0] == '+' or line[0] == '-')
    {
        sign = line[0];
        line.erase(0, 1);
    }

    size_t sz{0};
    for (; sz < line.size() and (std::isdigit(line[sz]) or line[sz] == '/');
         ++sz)
        ;

    /* parse the coefficient without the sign */
    std::string const tmp{line, 0, sz};
    fmpq_class cf = (sz == 0) ? fmpq_class(1, 1) : fmpq_class(tmp, 10);

    if (sign == '-')
        cf *= -1;

    cf.canonicalize();

    *pos = sz + (sign != 0);
    return cf;
}

}  // namespace

void generators_data::read_num_vars(std::FILE* infile)
{
    if (not getline(&m_line_ptr, &m_size, infile))
        throw std::runtime_error("Input file is empty.");

    std::string line{m_line_ptr};

    remove_spaces(line);
    remove_trailing_comma(line);

    /* the number of variables is the number of commas + 1 */
    num_vars = static_cast<uint32_t>(
        std::count(std::cbegin(line), std::cend(line), ',') + 1);
}

void generators_data::read_characteristic(std::FILE* infile)
{
    if (not getline(&m_line_ptr, &m_size, infile))
        throw std::runtime_error("No line containing field characteristic.");

    std::string line{m_line_ptr, m_size};

    uint64_t fc;
    try
    {
        fc = std::stoull(line);
    }
    catch (...)
    {
        throw std::runtime_error("Field characteristic not valid.");
    }

    if (fc < 0 or fc > std::numeric_limits<int32_t>::max())  // signed
    {
        throw std::runtime_error(
            "Field characteristic must be >= 0 and < 2^31.");
    }

    if (fc == 0)
        return;

    if (not n_is_prime(fc))
    {
        throw std::runtime_error("Field characteristic is not a prime number.");
    }

    field_char = static_cast<uint32_t>(fc);
}

void generators_data::read_num_generators(std::FILE* infile)
{
    num_gens = 0;

    while (getline(&m_line_ptr, &m_size, infile) != -1)
    {
        std::string line{m_line_ptr};

        /* check if there are empty lines in the input file */
        remove_spaces(line);
        remove_trailing_comma(line);

        if (not line.empty())
            ++num_gens;
    }
}

void generators_data::read_variable_names(std::FILE* infile)
{
    getline(&m_line_ptr, &m_size, infile);
    std::string line{m_line_ptr};

    remove_spaces(line);
    remove_trailing_comma(line);

    std::istringstream iline{line};

    /* variable names can't begin with one of this characters */
    static std::string const restricted_chars{"+-*^"};

    std::string var_name;
    while (std::getline(iline, var_name, ','))
    {
        if (var_name.empty())
        {
            throw std::runtime_error("Empty variable name.");
        }

        if (std::isdigit(var_name[0])
            or restricted_chars.find(var_name[0]) != std::string::npos)
        {
            throw std::runtime_error("Invalid variable name: " + var_name
                                     + ".");
        }

        if (std::find(std::cbegin(var_names), std::cend(var_names), var_name)
            != std::cend(var_names))
        {
            throw std::runtime_error("Duplicated variable name: " + var_name
                                     + ".");
        }

        var_names.emplace_back(var_name);
    }
}

void generators_data::read_exponent(std::string term)
{
    constexpr auto const npos = std::string::npos;

    static std::vector<int32_t> exp(num_vars);

    for (size_t k = 0; k < num_vars; ++k)
    {
        int32_t sexpn = 0; /* by default if varname is not yet found */

        /* use sorted vnames by length to avoid colisions with substrings */
        std::string const& vname = sorted_var_names[k];

        /* we allow repeated variables in the same monomial */
        for (size_t var = term.find(vname); var != npos; var = term.find(vname))
        {
            size_t next_var = var + vname.size();

            int32_t expn = 1; /* by default if varnmae is found */
            /* if there follows an exp symbol "^" */
            if (term[next_var] == '^')
            {
                size_t exp_size;
                expn = std::stoi(term.substr(++next_var), &exp_size, 10);

                constexpr int32_t const max_exponent =
                    1 << (CHAR_BIT * sizeof(exponent_type) - 1);

                if (expn < 0 or expn + sexpn >= max_exponent)
                {
                    throw std::out_of_range{"Exponent too large"};
                }

                next_var += exp_size;
            } /* else, the exponent is 1 by default and next_var is good to
               * continue */

            sexpn += expn;
            /* since variable names cannot start with digits we can or cannot
             * have the character '*' */
            if (term[next_var] == '*')
                ++next_var;

            /* erase monomial and exponent just parsed */
            term.erase(var, next_var - var);
        }

        /* save the total sum of the k-th exponent in the correct spot */
        exp[argsort_vnames[k]] = sexpn;
    }

    /* chars. left in 'term' after parsing monomials must be invalid input */
    if (not term.empty())
    {
        throw std::invalid_argument{"Error parsing monomial."};
    }

    exps.insert(std::end(exps), std::cbegin(exp), std::cend(exp));
}

void generators_data::read_generator_line(std::string line,
                                          size_t const line_num)
{
    size_t num_terms{0};

    for (auto it = std::begin(line); it != std::end(line);)
    {
        std::string term{read_term(line, it)};

        /* +-, -+, *x, 101*, x*y* are bad syntax */
        if (term == "+" or term == "-" or term.front() == '*'
            or term.back() == '*')
        {
            throw std::runtime_error("Error parsing generator number "
                                     + std::to_string(line_num) + ".");
        }

        size_t pos{0};
        fmpq_class cf;
        try
        {
            cf = term_to_mpq(term, &pos);
        }
        catch (std::invalid_argument&)
        {
            throw std::runtime_error("Error parsing generator number "
                                     + std::to_string(line_num) + ".");
        }

        /* denominators cannot be divisible by field characteristic */
        if (field_char != 0 and cf.get_den() % field_char == 0)
        {
            throw std::runtime_error("Division by zero in generator number "
                                     + std::to_string(line_num) + ".");
        }

        coeffs.emplace_back(std::move(cf));
        ++num_terms;

        /* skip the * in a monomial like 2*xy */
        if (term[pos] == '*')
            ++pos;

        try
        {
            generators_data::read_exponent(term.substr(pos));
        }
        catch (std::invalid_argument&)
        {
            throw std::runtime_error("Error parsing generator number "
                                     + std::to_string(line_num) + ".");
        }
        catch (std::out_of_range&)
        {
            constexpr int32_t const max_exp_bound =
                1 << (CHAR_BIT * sizeof(exponent_type) - 1);

            throw std::runtime_error(
                "Monomial exponent in generator " + std::to_string(line_num)
                + " must be >= 0 or < " + std::to_string(max_exp_bound) + ".");
        }
    }

    lens.emplace_back(num_terms);
}

void generators_data::read_generators(std::FILE* infile)
{
    /* ignore line with field characteristic */
    getline(&m_line_ptr, &m_size, infile);

    for (size_t lnum = 1; getline(&m_line_ptr, &m_size, infile) != -1; ++lnum)
    {
        std::string line{m_line_ptr};

        remove_spaces(line);
        remove_trailing_comma(line);

        /* skip empty lines in the input file */
        if (line.empty())
            continue;

        read_generator_line(line, lnum);
    }
}

void generators_data::read(std::FILE* infile)
{
    read_num_vars(infile);
    read_characteristic(infile);
    read_num_generators(infile);

    fseek(infile, 0, SEEK_SET);

    read_variable_names(infile);

    /* argsort variable names by descending length */
    argsort_vnames = sort_permutation(
        var_names, [](auto& v1, auto& v2) { return v1.size() > v2.size(); });
    /* use index permutation to get the actual sorted strings */
    sorted_var_names = apply_permutation(var_names, argsort_vnames);

    read_generators(infile);

    /* reduce coefficients modulo the field characteristic */
    if (field_char > 0)
    {
        coeffs_modp.resize(coeffs.size());

        std::ranges::transform(
            coeffs, std::begin(coeffs_modp), [this](fmpq_class const& c) {
                uint64_t const a =
                    fmpz_fdiv_ui(c.get_num().get_fmpz_t(), field_char);
                uint64_t const b =
                    fmpz_fdiv_ui(c.get_den().get_fmpz_t(), field_char);

                return n_mulmod2(a, n_invmod(b, field_char), field_char);
            });
    }

    free(m_line_ptr);
    m_line_ptr = nullptr;
    m_size     = 0ULL;
}

template <bool Star>
bool generators_data::write_monomial(std::FILE* outfile, size_t offset) const
{
    bool flag{true};

    for (size_t k = 0; k < num_vars; ++k)
    {
        auto const exp = exps[offset * num_vars + k];

        if (exp > 0)
        {
            if (Star or not flag)
                fmt::print(outfile, "*{}", var_names[k]);
            else
                fmt::print(outfile, "{}", var_names[k]);

            if (exp != 1)
                fmt::print(outfile, "^{}", exp);

            flag = false;
        }
    }

    return flag;
}

template <class CoefficientType>
void generators_data::write_generators(
    std::FILE* outfile,
    std::vector<CoefficientType> const& coeff) const
{
    for (size_t i = 0, offset = 0; i < num_gens; ++i)
    {
        for (size_t j = 0; j < lens[i]; ++j)
        {
            auto const& cf = coeff[offset + j];

            if (cf == 0)
                continue;

            if (cf != 1 and -cf != 1)
            {
                if constexpr (std::is_same_v<CoefficientType, fmpq_class>)
                {
                    flint_fprintf(outfile, "%{fmpz}",
                                  cf.get_num().get_fmpz_t());

                    if (cf.get_den() != 1)
                        flint_fprintf(outfile, "/%{fmpz}",
                                      cf.get_den().get_fmpz_t());
                }
                else
                    fmt::print(outfile, "{}", cf);

                write_monomial<true>(outfile, offset + j);
            }
            else
            {
                if (cf < 0)
                    fmt::print(outfile, "-");

                bool const flag = write_monomial<false>(outfile, offset + j);
                /* if no coefficient has been printed and the monomial is x^0
                 * print 1 */
                if (flag == 1)
                    fmt::print(outfile, "1");
            }

            if (j != lens[i] - 1 and coeff[offset + j + 1] > 0)
                fmt::print(outfile, "+");
        }

        if (i != num_gens - 1)
            fmt::print(outfile, ",");

        fmt::print(outfile, "\n");

        offset += lens[i];
    }
}

void generators_data::write(std::FILE* outfile) const
{
    for (auto const& vname : var_names)
        fmt::print(outfile, "{},", vname);

    fmt::print(outfile, "\n{}\n", field_char);

    if (field_char > 0)
    {
        write_generators(outfile, coeffs_modp);
    }
    else
    {
        write_generators(outfile, coeffs);
    }
}

}  // namespace gamba
