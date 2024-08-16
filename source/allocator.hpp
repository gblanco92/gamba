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

#include <cstddef>
#include <cstdlib>  // posix_memalign
#include <limits>
#include <new>
#include <type_traits>

#include <sys/mman.h>  // madvise

#include "config.hpp"

namespace gamba
{

/* Borrowed from: https://stackoverflow.com/questions/60169819/modern-approach-
 * to-making-stdvector-allocate-aligned-memory */

/*
 * Returns aligned pointers when allocations are requested. Default alignment
 * is 64B = 512b, sufficient for AVX-512 and most cache line sizes.
 *
 * @tparam AlignmentInBytes must be a positive power of 2.
 */
template <class ElementType, size_t const AlignmentInBytes = CACHE_LINE_SIZE>
class aligned_allocator
{
    static_assert(
        AlignmentInBytes >= alignof(ElementType),
        "Beware that types like int have minimum alignment requirements "
        "or access will result in crashes.");

public:
    using value_type = ElementType;

    static constexpr std::align_val_t const alignment{AlignmentInBytes};

    static constexpr std::size_t const huge_page_size{1UL << 21};  // 2 MiB

    /**
     * This is only necessary because aligned_allocator has a second template
     * argument for the alignment that will make the default
     * std::allocator_traits implementation fail during compilation.
     * @see https://stackoverflow.com/a/48062758/2191065
     */
    template <class OtherElementType>
    struct rebind
    {
        using other = aligned_allocator<OtherElementType, AlignmentInBytes>;
    };

    constexpr aligned_allocator() noexcept = default;

    constexpr aligned_allocator(aligned_allocator const&) noexcept = default;

    template <typename U>
    constexpr explicit aligned_allocator(
        aligned_allocator<U, AlignmentInBytes> const& /*unused*/) noexcept
    {}

    [[nodiscard]] ElementType* allocate(size_t num)
    {
        if (num > std::numeric_limits<size_t>::max() / sizeof(ElementType))
        {
            throw std::bad_array_new_length{};
        }

        void* ptr{nullptr};

        size_t const num_bytes = num * sizeof(ElementType);

        int error{0};

        if constexpr (AlignmentInBytes == huge_page_size)
        {
            error = posix_memalign(&ptr, huge_page_size, num_bytes);

#ifdef __linux__
            [[maybe_unused]] int const madvise_error =
                madvise(ptr, num_bytes, MADV_HUGEPAGE);

            // user's system may not enable kernel CONFIG_TRANSPARENT_HUGEPAGE
#    ifndef GAMBA_RELEASE
            if (madvise_error)
            {
                throw std::bad_alloc{};
            }
#    endif
#endif
        }
        else
        {
            ptr = ::operator new[](num_bytes, alignment);
        }

        if (error or ptr == nullptr)
        {
            throw std::bad_alloc{};
        }

        return reinterpret_cast<ElementType*>(ptr);
    }

    void deallocate(ElementType* ptr, [[maybe_unused]] size_t num_bytes = 0)
    {
        if constexpr (AlignmentInBytes == huge_page_size)
        {
            free(ptr);
        }
        else
        {
            /* According to the C++20 draft n4868 § 17.6.3.3, the delete
             * operator must be called with the same alignment argument as the
             * new expression. The size argument can be omitted but if present
             * must also be equal to the one used in new. */
            ::operator delete[](ptr, alignment);
        }
    }
}; /* end class aligned_allocator */

template <class T1, size_t N1, class T2, size_t N2>
constexpr bool operator==(aligned_allocator<T1, N1> const& /*unused*/,
                          aligned_allocator<T2, N2> const& /*unused*/) noexcept
{
    return (std::is_same_v<T1, T2> and N1 == N2);
}

template <class T1, size_t N1, class T2, size_t N2>
constexpr bool operator!=(aligned_allocator<T1, N1> const& lhs,
                          aligned_allocator<T2, N2> const& rhs) noexcept
{
    return not(lhs == rhs);
}

}  // namespace gamba
