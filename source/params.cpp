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

#include "params.hpp"

namespace gamba
{

std::map<std::string, order_options> gamba_params::mon_order_map = {
    {  "grevlex",   order_options::grevlex},
    {    "lexic",     order_options::lexic},
    {"blockelim", order_options::blockelim}
};

void check_gamba_parameters(gamba_params& params,
                            generators_data const& input_data)
{
    /* zero means taking all the spairs with same degree */
    if (params.max_spairs == 0UL)
    {
        params.max_spairs = std::numeric_limits<ssize_t>::max();
    }

    /* validate num_elim_vars & set monomial order */
    if (params.num_elim_vars != 0)
    {
        if (params.num_elim_vars >= input_data.num_vars)
        {
            throw std::runtime_error("Elimination block size must be "
                                     "smaller than the number of variables.");
        }

        params.mon_order = order_options::blockelim;
    }

    /* select what information to display depending on verbosity level */
    params.basis_info   = params.verbose > 0;
    params.rounds_info  = params.verbose > 1;
    params.timings_info = params.verbose > 0;
    params.stats_info   = params.verbose > 0;
}

}  // namespace gamba
