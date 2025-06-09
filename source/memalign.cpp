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

#include "memalign.hpp"

#include <sys/mman.h>

#include "config.hpp"
#include "utils.hpp"

namespace gamba
{

void* memalign_alloc(size_t const size) noexcept
{
    static constexpr size_t const alignment = HUGE_PAGE_SIZE;
    /* alignment is a power of 2 */
    static_assert((alignment & (alignment - 1)) == 0);
    /* don't do stupid allocations; mmap will fail anyway */
    assert(size != 0UL);

    /* for all requested memory to reside in huge pages we must work with whole
     * huge pages chunks (otherwise the remainder will not be in a huge page) */
    size_t const rounded_size = round_up(size, alignment);
    /* mmap returns page-aligned memory so an extra huge page chunk is needed
     * to fit all the (rounded) requested memory in huge pages */
    size_t const mmap_size = rounded_size + alignment;

    /* we want transparent huge pages; MAP_HUGETLB uses permanent huge pages */
    void* const block = ::mmap(nullptr, mmap_size, PROT_READ | PROT_WRITE,
                               MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    if (block == MAP_FAILED)
    {
        return nullptr;
    }

    auto const block_addr = reinterpret_cast<std::uintptr_t>(block);
    /* get pointer by aligning memory block upwards */
    std::uintptr_t const ptr_addr = (block_addr + (alignment - 1)) & -alignment;

    /* check that pointer alignment is always correct */
    assert(ptr_addr % alignment == 0);

    size_t const offset   = block_addr % alignment;
    size_t const head_gap = (offset != 0 ? alignment - offset : 0UL);
    size_t const tail_gap = alignment - head_gap;
    assert(head_gap + rounded_size + tail_gap == mmap_size);

    std::uintptr_t const block_end_addr = block_addr + head_gap + rounded_size;
    auto* const block_end = reinterpret_cast<void*>(block_end_addr);  // NOLINT

    /* munmap head/tail gaps; doing this is the only way mremap respects the
     * huge pages (and the newly allocated memory is also huge pages) */
    if (head_gap != 0)
        ::munmap(block, head_gap);
    if (tail_gap != 0)
        ::munmap(block_end, tail_gap);

    /* at this point the mmap'ed size is the requested (rounded) size */
    auto* const ptr = reinterpret_cast<void*>(ptr_addr);  // NOLINT

    /* only Linux and some BSD kernels have transparent huge pages support */
#ifdef MADV_HUGEPAGE
    /* use madvise to enable transparent huge pages */
    [[maybe_unused]] int const madvise_error =
        ::madvise(ptr, rounded_size, MADV_HUGEPAGE);

    assert(not madvise_error);
    /* do not check for madvise error as user's system may not enable
     * transparent huge pages by default */
#endif

    return ptr;
}

void memalign_free(void* const ptr, size_t const size) noexcept
{
    if (ptr == nullptr)
        return;

    /* stupid allocations are not possible */
    assert(size != 0UL);

    /* memalign functions internally allocate multiples of HUGE_PAGE_SIZE */
    [[maybe_unused]] int const error =
        ::munmap(ptr, round_up(size, HUGE_PAGE_SIZE));

    assert(not error);
}

}  // namespace gamba
