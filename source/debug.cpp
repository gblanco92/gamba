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

#include "debug.hpp"

namespace gamba
{

void print_compiler_version(std::ostream& out)
{
    out << "Compiler version: " << cxx << "-" << cxx_ver_major << "."
        << cxx_ver_minor << "." << cxx_ver_patch << std::endl;
}

void print_libgmp_version(std::ostream& out)
{
    out << "libGMP version: " << __GNU_MP_VERSION << "."
        << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL
        << std::endl;
}

int32_t line_inserter::operator()(std::streambuf& sbuf, int32_t const c)
{
    int32_t const result = eof();
    if (c != eof() and m_new_line)
        std::ostream(&sbuf) << "## DEBUG ## ";

    m_new_line = (c == '\n');

    assert(0 <= c and c <= std::numeric_limits<char>::max());
    return result != eof() ? result : sbuf.sputc(static_cast<char>(c));
}

}  // namespace gamba
