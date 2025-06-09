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

#include "utils.hpp"

namespace gamba
{

struct divmask_type
{
    using mask_type = uint64_t;

    mask_type mask;

    /* if the monomial x_a with divmask a divides the monomial x_b with divmask
     * b, then a <= b, hence !(a <= b) implies x_a doesn't divide x_b */
    FORCE_INLINE bool operator<=(divmask_type other) const
    {
        return (this->mask & ~other.mask) == 0;
    }
};

static_assert(std::is_trivial_v<divmask_type>);

}  // namespace gamba
