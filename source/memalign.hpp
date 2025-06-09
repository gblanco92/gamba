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

/*
 * The function 'memalign_alloc' returns a pointer to a memory region of at
 * least 'size' bytes that is hugepage-aligned and already marked as hugepage,
 * if the underlying OS supports transparent huge pages via 'madvise'. Any
 * allocation performed via 'memalign_alloc' must be deallocated using
 * 'memaling_free'. 'memalign_alloc' and 'memalign_realloc' return 'nullptr' if
 * an error occurs. In constrast to 'posix_memalign' all memory is allocated
 * directly using mmap and all memory (even if it is smaller than a single huge
 * page) will be placed inside a (transparent) huge page.
 */

[[nodiscard]] void* memalign_alloc(size_t const size) noexcept;

void memalign_free(void* ptr, size_t const size) noexcept;

}  // namespace gamba
