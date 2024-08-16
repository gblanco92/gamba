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

#include "io.hpp"

namespace gamba
{

enum class order_options
{
    grevlex,
    lexic,
    blockelim
};

struct gamba_params
{
    bool all_spairs{false};

    ssize_t max_spairs{2'000L};

    order_options mon_order{order_options::grevlex};

    static std::map<std::string, order_options> mon_order_map;

    size_t num_elim_vars{0UL};

    size_t num_threads{1UL};

    bool no_reduce{false};

    size_t seed{967'557'673UL};

    size_t verbose{2};

    bool basis_info{true};

    bool rounds_info{true};

    bool stats_info{true};

    bool timings_info{true};
};

void check_gamba_parameters(gamba_params& params,
                            generators_data const& input_data);

}  // namespace gamba
