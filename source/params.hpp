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

#pragma once

namespace gamba
{

struct params
{
    enum class order  // NOLINT
    {
        grevlex,
        deglex,
        lexic,
        blockelim,
        grevlexw,
    };

    static void sanitize_input(size_t const num_vars);

    /* GamBa parameters */

    static bool all_spairs;

    static ssize_t max_spairs;

    static std::string mon_order_str;

    static order mon_order;

    static size_t num_elim_vars;

    static std::vector<uint32_t> weights;

    static bool lead_mons;

    static size_t num_threads;

    static bool no_reduce;

    static size_t seed;

    static ssize_t verbose;
};

}  // namespace gamba
