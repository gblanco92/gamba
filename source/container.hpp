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

#include "allocator.hpp"
#include "config.hpp"
#include "flat_set.hpp"

namespace gamba
{

template <class MonomialType>
using monomial_set = monomial_flat_set<MonomialType>;

template <class ElementType>
using aligned_vector =
    std::vector<ElementType, aligned_allocator<ElementType, CACHE_LINE_SIZE>>;

constexpr size_t const HUGE_PAGE_SIZE = aligned_allocator<int>::huge_page_size;

template <class ElementType>
using hugepage_vector =
    std::vector<ElementType, aligned_allocator<ElementType, HUGE_PAGE_SIZE>>;

}  // namespace gamba
