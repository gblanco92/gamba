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

#include <stdexcept>

#include "params.hpp"

#include "utils.hpp"

namespace gamba
{

bool params::all_spairs{false};

ssize_t params::max_spairs{2'000L};

std::string params::mon_order_str{"grevlex"};

params::order params::mon_order{params::order::grevlex};

size_t params::num_elim_vars{0UL};

std::vector<uint32_t> params::weights{};

bool params::lead_mons{false};

size_t params::num_threads{1UL};

bool params::no_reduce{false};

GAMBA_RELEASE(size_t params::seed = std::random_device{}();)

GAMBA_DEVELOP(size_t params::seed{967'557'673UL};)

ssize_t params::verbose{2};

void params::sanitize_input(size_t const num_vars)
{
    /* zero means taking all the spairs with same degree */
    if (params::max_spairs == 0UL)
    {
        params::max_spairs = std::numeric_limits<ssize_t>::max();
    }

    /* set monomial order from input string */
    if (params::mon_order_str == "grevlex")
        params::mon_order = params::order::grevlex;
    else if (params::mon_order_str == "deglex")
        params::mon_order = params::order::deglex;
    else if (params::mon_order_str == "lexic")
        params::mon_order = params::order::lexic;
    else if (params::mon_order_str == "blockelim")
        params::mon_order = params::order::blockelim;
    else if (params::mon_order_str == "grevlexw")
        params::mon_order = params::order::grevlexw;
    else
        throw std::runtime_error("Invalid monomial order.");

    /* if num_elim_vars is different from zero 'blockelim' must be used */
    if (params::num_elim_vars != 0)
    {
        if (params::mon_order != params::order::blockelim)
            throw std::runtime_error("Elimination of variables requires "
                                     "'blockelim' monomial order.");
    }

    /* if 'blockelim' is used check num_elim_vars is not too big */
    if (params::mon_order == params::order::blockelim)
    {
        if (params::num_elim_vars <= 0)
            throw std::runtime_error(
                "Elimination block size must be non-zero.");

        if (params::num_elim_vars >= num_vars)
            throw std::runtime_error("Elimination block size must be "
                                     "smaller than the number of variables.");
    }

    /* if weights has been set, 'grevlexw' must be used */
    if (not params::weights.empty())
    {
        if (params::mon_order != params::order::grevlexw)
            throw std::runtime_error(
                "Monomial weights requires 'grevlexw' monomial order.");
    }

    /* if 'grevlex' is used check that we have the right number of weights */
    if (params::mon_order == params::order::grevlexw)
    {
        if (params::weights.size() != num_vars)
            throw std::runtime_error("Monomial weights length must be equal to "
                                     "the number of variables.");
    }
}

}  // namespace gamba
