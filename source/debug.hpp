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

#include <cassert>
#include <iostream>

#include <gmpxx.h>

#ifdef DEBUG
#    define GAMBA_DEBUG(code) code
#else
#    define GAMBA_DEBUG(code)
#endif

#ifdef DEBUG_STDLIB
#    if defined(__GLIBCXX__) && DEBUG_STDLIB == 0
#        define _GLIBXCXX_ASSERTIONS
#    elif defined(__GLIBCXX__) && DEBUG_STDLIB >= 1
#        define _GLIBCXX_DEBUG
#    elif defined(_LIBCPP_VERSION) && DEBUG_STDLIB == 0
#        define _LIBCPP_DEBUG 0
#    elif defined(_LIBCPP_VERSION) && DEBUG_STDLIB >= 1
#        define _LIBCPP_DEBUG 1
#    endif
#endif

namespace gamba
{

constexpr std::string_view const cxx =
#ifdef __clang__
    "clang++";
#else
    "g++";
#endif

constexpr int32_t const cxx_ver_major =
#ifdef __clang__
    __clang_major__;
#else
    __GNUC__;
#endif

constexpr int32_t const cxx_ver_minor =
#ifdef __clang__
    __clang_minor__;
#else
    __GNUC_MINOR__;
#endif

constexpr int32_t const cxx_ver_patch =
#ifdef __clang__
    __clang_patchlevel__;
#else
    __GNUC_PATCHLEVEL__;
#endif

void print_compiler_version(std::ostream& out);

void print_libgmp_version(std::ostream& out);

/* Borrowed from: https://stackoverflow.com/questions/23688447/inserting-text-
 * before-each-line-using-stdostream */
template <class Inserter>
class filtering_streambuf : public std::streambuf
{
public:
    filtering_streambuf(std::streambuf* sbuf,
                        Inserter inserter,
                        bool delete_when_finished) :
            m_sbuf{sbuf},
            m_inserter{inserter},
            m_delete_when_finished{delete_when_finished}
    {}

    explicit filtering_streambuf(std::streambuf* sbuf,
                                 bool delete_when_finished = false) :
            m_sbuf{sbuf},
            m_inserter{},
            m_delete_when_finished{delete_when_finished}
    {}

    filtering_streambuf(filtering_streambuf const&)            = delete;
    filtering_streambuf& operator=(filtering_streambuf const&) = delete;

    ~filtering_streambuf() override
    {
        sync();

        if (m_delete_when_finished)
            delete m_sbuf;
    }

    int_type overflow(int_type c) override
    {
        if (c == traits_type::eof())
            sync();
        else if (m_sbuf != nullptr)
            m_inserter(*m_sbuf, c);

        return c;
    }

    int32_t sync() override { return std::streambuf::sync(); }

private:
    std::streambuf* m_sbuf;

    Inserter m_inserter;

    bool m_delete_when_finished;
};

class line_inserter
{
public:
    line_inserter() = default;

    int32_t operator()(std::streambuf& sbuf, int32_t const c);

private:
    static constexpr int32_t eof() { return std::char_traits<char>::eof(); }

    bool m_new_line{true};
};

static filtering_streambuf<line_inserter> linebuf(std::cout.rdbuf());  // NOLINT
static std::ostream debug(&linebuf);

}  // namespace gamba
